'''
======================================
Expression Dotplot Colored by Survival
======================================

Adapted from vivarium-cell.
'''
from matplotlib import pyplot as plt
import numpy as np

from vivarium.core.experiment import get_in
from vivarium.core.emitter import path_timeseries_from_data

from src.investigate_utils import filter_raw_data_by_time


PATH_TO_AGENTS = ('agents',)
PATH_TO_DEAD = ('boundary', 'dead')
LIVE_COLOR = 'green'
DEAD_COLOR = 'black'
ALPHA = 0.5


def _get_final_live_agents(data, time_range=(0, 1)):
    data = filter_raw_data_by_time(data, time_range)
    max_time = max(data.keys())
    agents = []
    agents_data = get_in(data[max_time], PATH_TO_AGENTS)
    for agent, agent_data in agents_data.items():
        dead = get_in(agent_data, PATH_TO_DEAD)
        if not dead:
            agents.append(agent)
    return agents


def plot_expression_survival(
    data, path_to_x_variable, path_to_y_variable, xlabel, ylabel,
    boundary_x, boundary_y, boundary_error, boundary_color='black',
    scaling=1, time_range=(0, 1), label_agents=False,
    plot_agents=tuple(), fontsize=36,
    dead_trace_agents=tuple(), agents_for_phylogeny_trace=tuple(),
):
    '''Create Expression Scatterplot Colored by Survival

    Plot one dot for each cell along an axis to indicate that cell's
    final expression level for a specified protein. The dot color
    reflects whether the cell survived long enough to divide.

    Parameters:
        data (dict): The raw data emitted from the simulation.
        path_to_x_variable (tuple): Path from the agent root to the
            variable that holds the protein's expression level. We do
            not adjust for cell volume, so this should be a
            concentration. This protein will be plotted on the x axis.
        path_to_y_variable (tuple): Path from the agent root to the
            variable that holds the protein's expression level. We do
            not adjust for cell volume, so this should be a
            concentration. This protein will be plotted on the y axis.
        xlabel (str): Label for x-axis.
        xlabel (str): Label for y-axis.
        boundary_x (list(float)): X-values of the boundary identified
            numerically from the antibiotic model.
        boundary_y (list(float)): Y-values of the boundary identified
            numerically from the antibiotic model.
        boundary_error (list(float)): Precision of the Y-value
            predictions. This is the distance along the y axis between
            the known dead point and known live point closest to the
            prediction. The predicted y-value will be in the middle of
            this range, so an error band can be drawn of width
            boundary_error centered at boundary_y.
        boundary_color (str): Color of boundary.
        scaling (str): Coefficient to multiply all data by. This is
            intended to be used for changing the units plotted.
        time_range (tuple): Tuple of two :py:class:`float`s that are
            fractions of the total simulated time period. These
            fractions indicate the start and end points (inclusive) of
            the time range to consider.
        label_agents (bool): Whether to label each point with the agent
            ID.
        dead_trace_agents (Iterable): The agent IDs of the agents to
            plot traces for. By default, no traces are shown.
        plot_agents(Iterable): The agent IDs of the agents to plot. By
            default, all agents are plotted.
        agents_for_phylogeny_trace (Iterable): Agent IDs for the agents
            whose phylogenies will be traced.
        fontsize (float): Text size for entire figure.

    Returns:
        plt.Figure: The finished figure.
    '''
    live_finals_x, dead_finals_x = calc_live_and_dead_finals(
        data, path_to_x_variable, time_range, plot_agents)
    live_finals_y, dead_finals_y = calc_live_and_dead_finals(
        data, path_to_y_variable, time_range, plot_agents)
    if label_agents:
        fig, ax = plt.subplots(figsize=(50, 50))
        # Always trace all agents when labeling
        dead_trace_agents = (
            list(dead_finals_x.keys()))
        agents_for_phylogeny_trace = _get_final_live_agents(
            data, time_range)
    else:
        fig, ax = plt.subplots(figsize=(6.4, 7))

    ax.scatter(
        np.array(list(live_finals_x.values())) * scaling,
        np.array(list(live_finals_y.values())) * scaling,
        label='Final Location at Division (Cell Survives)',
        color=LIVE_COLOR, alpha=ALPHA,
    )
    ax.scatter(
        np.array(list(dead_finals_x.values())) * scaling,
        np.array(list(dead_finals_y.values())) * scaling,
        label='Final Location at Death', color=DEAD_COLOR, alpha=ALPHA,
    )
    if label_agents:
        for agent in live_finals_x:
            x = live_finals_x[agent] * scaling
            y = live_finals_y[agent] * scaling
            ax.annotate(agent, (x, y), size=0.1)
        for agent in dead_finals_x:
            x = dead_finals_x[agent] * scaling
            y = dead_finals_y[agent] * scaling
            ax.annotate(agent, (x, y), size=0.1)
    plot_expression_survival_traces(ax, data, path_to_x_variable,
            path_to_y_variable, scaling, time_range, dead_trace_agents,
            DEAD_COLOR, agents_for_phylogeny_trace, LIVE_COLOR)
    finals = list(live_finals_x.values()) + list(
        dead_finals_x.values())
    boundary_x_arr = np.array(boundary_x)
    boundary_y_arr = np.array(boundary_y)
    mask = (
        (min(finals) <= boundary_x_arr)
        & (boundary_x_arr <= max(finals)))
    true_indices = np.where(mask)[0]
    if len(true_indices) > 0:
        min_true_idx = min(true_indices)
        max_true_idx = max(true_indices)
        # Make sure boundary spans entire figure
        if min_true_idx > 0:
            mask[min_true_idx - 1] = True
        if max_true_idx < len(mask) - 1:
            mask[max_true_idx + 1] = True
    boundary_x_arr = boundary_x_arr[mask]
    boundary_y_arr = boundary_y_arr[mask]
    ax.plot(
        boundary_x_arr * scaling, boundary_y_arr * scaling,
        c=boundary_color,
        label='Boundary Predicted by Antibiotics Model')
    ax.fill_between(
        boundary_x_arr * scaling,
        boundary_y_arr * scaling - boundary_error * scaling / 2,
        boundary_y_arr * scaling + boundary_error * scaling / 2,
        color=boundary_color, alpha=0.2)
    ax.legend(  # type: ignore
        bbox_to_anchor=(0.5, 1.05), loc='lower center',
        prop={'size': fontsize})
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    for spine_name in ('top', 'right'):
        ax.spines[spine_name].set_visible(False)
    fig.tight_layout()
    return fig


def plot_expression_survival_traces(
    ax, data, path_to_x_variable, path_to_y_variable, scaling=1,
    time_range=(0, 1), dead_agents=tuple(), dead_trace_color='black',
    agents_for_phylogeny_trace=tuple(), phylogeny_trace_color='green',
):
    '''Create Expression Traces Colored by Survival

    Plot a trace of expression values for each cell. The color of each
    dot in the trace reflects whether the cell was alive at that point.

    Parameters:
        ax (Axes): Axes to plot on.
        data (dict): The raw data emitted from the simulation.
        path_to_x_variable (tuple): Path from the agent root to the
            variable that holds the protein's expression level. We do
            not adjust for cell volume, so this should be a
            concentration. This protein will be plotted on the x axis.
        path_to_y_variable (tuple): Path from the agent root to the
            variable that holds the protein's expression level. We do
            not adjust for cell volume, so this should be a
            concentration. This protein will be plotted on the y axis.
        scaling (str): Coefficient to multiply all data by. This is
            intended to be used for changing the units plotted.
        time_range (tuple): Tuple of two :py:class:`float`s that are
            fractions of the total simulated time period. These
            fractions indicate the start and end points (inclusive) of
            the time range to consider.
        dead_agents (Iterable): The agent IDs of the agents to plot
            traces for. These agents should die.
        dead_trace_color (str): Color of trace line for dead cells.
        agents_for_phylogeny_trace (Iterable): Agent IDs for the agents
            whose phylogenies will be traced.
        phylogeny_trace_color (str): Color of trace line for phylogeny.
    '''
    data = filter_raw_data_by_time(data, time_range)
    path_timeseries = path_timeseries_from_data(data)

    # Plot dead traces
    for i, agent in enumerate(dead_agents):
        x_timeseries = path_timeseries[
            PATH_TO_AGENTS + (agent,)
            + path_to_x_variable]
        y_timeseries = path_timeseries[
            PATH_TO_AGENTS + (agent,)
            + path_to_y_variable]
        ax.plot(
            np.array(x_timeseries) * scaling,
            np.array(y_timeseries) * scaling,
            color=dead_trace_color,
            linewidth=0.5,
            label='Agent path until death' if i == 0 else '',
        )

    # Plot phylogeny traces
    plotted_solid = False
    plotted_dashed = False
    for agent in agents_for_phylogeny_trace:
        last_end_point = tuple()
        for i in range(len(agent) + 1):
            ancestor = agent[:i]
            x_path = PATH_TO_AGENTS + (ancestor,) + path_to_x_variable
            y_path = PATH_TO_AGENTS + (ancestor,) + path_to_y_variable
            if x_path in path_timeseries:
                # else ancestor does not exist
                x_timeseries = path_timeseries[x_path]
                y_timeseries = path_timeseries[y_path]
                ax.plot(
                    np.array(x_timeseries) * scaling,
                    np.array(y_timeseries) * scaling,
                    color=phylogeny_trace_color,
                    linewidth=0.5,
                    label=(
                        'Agent path until division'
                        if not plotted_solid else ''
                    ),
                )
                if not plotted_solid:
                    # This is the first agent in the lineage
                    ax.scatter(
                        x_timeseries[0] * scaling,
                        y_timeseries[0] * scaling,
                        color=phylogeny_trace_color,
                        marker='s',
                        label='Lineage start',
                    )
                plotted_solid = True
                if last_end_point:
                    ax.plot(
                        [
                            last_end_point[0] * scaling,
                            x_timeseries[0] * scaling
                        ],
                        [
                            last_end_point[1] * scaling,
                            y_timeseries[0] * scaling
                        ],
                        color=phylogeny_trace_color,
                        linewidth=0.5,
                        linestyle='--',
                        label=(
                            'From final mother location to initial '
                            'daughter location'
                            if not plotted_dashed else ''
                        ),
                    )
                    plotted_dashed = True
                last_end_point = x_timeseries[-1], y_timeseries[-1]


def plot_expression_survival_dotplot(
    data, path_to_variable, xlabel, scaling=1, time_range=(0, 1),
    fontsize=36
):
    '''Create Expression Dotplot Colored by Survival

    Plot one dot for each cell along an axis to indicate that cell's
    final expression level for a specified protein. The dot color
    reflects whether the cell survived long enough to divide.

    Note that only the expression levels while the cell is alive are
    considered.

    Parameters:
        data (dict): The raw data emitted from the simulation.
        path_to_variable (tuple): Path from the agent root to the
            variable that holds the protein's expression level. We do
            not adjust for cell volume, so this should be a
            concentration.
        xlabel (str): Label for x-axis.
        scaling (str): Coefficient to multiply all data by. This is
            intended to be used for changing the units plotted.
        time_range (tuple): Tuple of two :py:class:`float`s that are
            fractions of the total simulated time period. These
            fractions indicate the start and end points (inclusive) of
            the time range to consider.
        fontsize (float): size of all text in figure.

    Returns:
        plt.Figure: The finished figure.
    '''
    live_finals, dead_finals = calc_live_and_dead_finals(
        data, path_to_variable, time_range)
    fig, ax = plt.subplots(figsize=(6, 2))
    ax.scatter(
        np.array(list(live_finals.values())) * scaling,
        [0.1] * len(live_finals),
        label='Survive', color=LIVE_COLOR, alpha=ALPHA,
    )
    ax.scatter(
        np.array(list(dead_finals.values())) * scaling,
        [0.1] * len(dead_finals),
        label='Die', color=DEAD_COLOR, alpha=ALPHA,
    )
    ax.legend(prop={'size': fontsize})
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylim([0, 1.25])
    ax.get_yaxis().set_visible(False)
    for spine_name in ('left', 'top', 'right'):
        ax.spines[spine_name].set_visible(False)
    ax.spines['bottom'].set_position('zero')
    fig.tight_layout()
    stats = {
        'live': list(live_finals.values()),
        'dead': list(dead_finals.values()),
    }
    return fig, stats


def calc_live_and_dead_finals(
        data, path_to_variable, time_range, agents=tuple()):
    values = {}
    die = set()

    end_time = max(data.keys())
    for time, time_data in data.items():
        if (time < time_range[0] * end_time
                or time > time_range[1] * end_time):
            continue
        agents_data = get_in(time_data, PATH_TO_AGENTS)
        for agent, agent_data in agents_data.items():
            if agents and agent not in agents:
                continue
            agent_values = values.setdefault(agent, {})
            value = get_in(agent_data, path_to_variable)
            if get_in(agent_data, PATH_TO_DEAD, False):
                die.add(agent)
            agent_values[time] = value

    live_finals = {}
    dead_finals = {}
    for agent, agent_values in values.items():
        if not agent_values:
            continue
        if agent in die:
            dead_finals[agent] = agent_values[max(agent_values.keys())]
        else:
            live_finals[agent] = agent_values[max(agent_values.keys())]
    return live_finals, dead_finals
