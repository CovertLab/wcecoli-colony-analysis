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
MARKERS = ('.', 'v', '^', 's', 'p', '*', '+', 'x', 'D')
ALPHA = 0.5


def plot_expression_survival(
    data, path_to_x_variable, path_to_y_variable, xlabel, ylabel,
    boundary_x, boundary_y, boundary_error, boundary_color='black',
    scaling=1, time_range=(0, 1), label_agents=False,
    trace_agents=tuple(),
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
        trace_agents (Iterable): The agent IDs of the agents to plot
            traces for.

    Returns:
        plt.Figure: The finished figure.
    '''
    live_finals_x, dead_finals_x = calc_live_and_dead_finals(
        data, path_to_x_variable, time_range)
    live_finals_y, dead_finals_y = calc_live_and_dead_finals(
        data, path_to_y_variable, time_range)
    fig, ax = plt.subplots()
    ax.scatter(
        np.array(list(live_finals_x.values())) * scaling,
        np.array(list(live_finals_y.values())) * scaling,
        label='Survive', color=LIVE_COLOR, alpha=ALPHA,
    )
    ax.scatter(
        np.array(list(dead_finals_x.values())) * scaling,
        np.array(list(dead_finals_y.values())) * scaling,
        label='Die', color=DEAD_COLOR, alpha=ALPHA,
    )
    if label_agents:
        for agent in live_finals_x:
            x = live_finals_x[agent] * scaling
            y = live_finals_y[agent] * scaling
            ax.annotate(agent, (x, y), size=1)
        for agent in dead_finals_x:
            x = dead_finals_x[agent] * scaling
            y = dead_finals_y[agent] * scaling
            ax.annotate(agent, (x, y), size=1)
    plot_expression_survival_traces(ax, data, path_to_x_variable,
            path_to_y_variable, scaling, time_range, trace_agents,
            LIVE_COLOR)
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
    ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    for spine_name in ('top', 'right'):
        ax.spines[spine_name].set_visible(False)
    fig.tight_layout()
    return fig


def plot_expression_survival_traces(
    ax, data, path_to_x_variable, path_to_y_variable, scaling=1,
    time_range=(0, 1), agents=tuple(), trace_color='black',
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
        agents (Iterable): The agent IDs of the agents to plot traces
            for.
        trace_color (str): Color of trace line.
    '''
    data = filter_raw_data_by_time(data, time_range)
    path_timeseries = path_timeseries_from_data(data)
    for agent in agents:
        x_timeseries = path_timeseries[
            PATH_TO_AGENTS + (agent,)
            + path_to_x_variable]
        y_timeseries = path_timeseries[
            PATH_TO_AGENTS + (agent,)
            + path_to_y_variable]
        ax.plot(
            np.array(x_timeseries) * scaling,
            np.array(y_timeseries) * scaling,
            color=trace_color,
            linewidth=0.5)


def plot_expression_survival_dotplot(
    data, path_to_variable, xlabel, scaling=1, time_range=(0, 1)
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
    ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylim([0, 1.25])
    ax.get_yaxis().set_visible(False)
    for spine_name in ('left', 'top', 'right'):
        ax.spines[spine_name].set_visible(False)
    ax.spines['bottom'].set_position('zero')
    fig.tight_layout()
    return fig


def calc_live_and_dead_finals(data, path_to_variable, time_range):
    values = {}
    die = set()

    end_time = max(data.keys())
    for time, time_data in data.items():
        if (time < time_range[0] * end_time
                or time > time_range[1] * end_time):
            continue
        agents_data = get_in(time_data, PATH_TO_AGENTS)
        for agent, agent_data in agents_data.items():
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
