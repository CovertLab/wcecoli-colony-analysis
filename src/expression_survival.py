'''
======================================
Expression Dotplot Colored by Survival
======================================

Adapted from vivarium-cell.
'''
from bisect import bisect_left, bisect_right

from matplotlib import pyplot as plt
import numpy as np

from vivarium.core.experiment import get_in
from vivarium.core.emitter import path_timeseries_from_data


PATH_TO_AGENTS = ('agents',)
PATH_TO_DEAD = ('boundary', 'dead')
LIVE_COLOR = 'green'
DEAD_COLOR = 'black'
MARKERS = ('.', 'v', '^', 's', 'p', '*', '+', 'x', 'D')
ALPHA = 0.5


def plot_expression_survival(
    data, path_to_x_variable, path_to_y_variable, xlabel, ylabel,
    boundary_x, boundary_y, boundary_error, scaling=1,
    time_range=(0, 1), label_agents=False,
):
    '''Create Expression Scatterplot Colored by Survival

    Plot one dot for each cell along an axis to indicate that cell's
    average expression level for a specified protein. The dot color
    reflects whether the cell survived long enough to divide.

    Note that only the expression levels while the cell is alive are
    considered in the average.

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
        scaling (str): Coefficient to multiply all data by. This is
            intended to be used for changing the units plotted.
        time_range (tuple): Tuple of two :py:class:`float`s that are
            fractions of the total simulated time period. These
            fractions indicate the start and end points (inclusive) of
            the time range to consider when calculating average
            expression level.
        label_agents (bool): Whether to label each point with the agent
            ID.

    Returns:
        plt.Figure: The finished figure.
    '''
    live_averages_x, dead_averages_x = calc_live_and_dead_averages(
        data, path_to_x_variable, time_range)
    live_averages_y, dead_averages_y = calc_live_and_dead_averages(
        data, path_to_y_variable, time_range)
    fig, ax = plt.subplots()
    ax.scatter(
        np.array(list(live_averages_x.values())) * scaling,
        np.array(list(live_averages_y.values())) * scaling,
        label='Survive', color=LIVE_COLOR, alpha=ALPHA,
    )
    ax.scatter(
        np.array(list(dead_averages_x.values())) * scaling,
        np.array(list(dead_averages_y.values())) * scaling,
        label='Die', color=DEAD_COLOR, alpha=ALPHA,
    )
    if label_agents:
        for agent in live_averages_x:
            x = live_averages_x[agent] * scaling
            y = live_averages_y[agent] * scaling
            ax.annotate(agent, (x, y), size=1)
        for agent in dead_averages_x:
            x = dead_averages_x[agent] * scaling
            y = dead_averages_y[agent] * scaling
            ax.annotate(agent, (x, y), size=1)
    averages = list(live_averages_x.values()) + list(
        dead_averages_x.values())
    boundary_x_arr = np.array(boundary_x)
    boundary_y_arr = np.array(boundary_y)
    mask = (
        (min(averages) <= boundary_x_arr)
        & (boundary_x_arr <= max(averages)))
    boundary_x_arr = boundary_x_arr[mask]
    boundary_y_arr = boundary_y_arr[mask]
    ax.plot(
        boundary_x_arr * scaling, boundary_y_arr * scaling, c='black',
        label='Boundary Predicted by Antibiotics Model')
    ax.fill_between(
        boundary_x_arr * scaling,
        boundary_y_arr * scaling - boundary_error * scaling / 2,
        boundary_y_arr * scaling + boundary_error * scaling / 2,
        color='black', alpha=0.2)
    ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    for spine_name in ('top', 'right'):
        ax.spines[spine_name].set_visible(False)
    fig.tight_layout()
    return fig


def plot_expression_survival_traces(
    data, path_to_x_variable, path_to_y_variable, xlabel, ylabel,
    scaling=1, time_range=(0, 1), agents=tuple(),
):
    '''Create Expression Traces Colored by Survival

    Plot a trace of expression values for each cell. The color of each
    dot in the trace reflects whether the cell was alive at that point.

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
        scaling (str): Coefficient to multiply all data by. This is
            intended to be used for changing the units plotted.
        time_range (tuple): Tuple of two :py:class:`float`s that are
            fractions of the total simulated time period. These
            fractions indicate the start and end points (inclusive) of
            the time range to consider when calculating average
            expression level.
        agents (Iterable): The agent IDs of the agents to plot traces
            for.

    Returns:
        plt.Figure: The finished figure.
    '''
    path_timeseries = path_timeseries_from_data(data)
    fig, ax = plt.subplots()
    times = path_timeseries['time']
    min_idx = bisect_left(times, time_range[0] * times[-1])
    max_idx = bisect_right(times, time_range[1] * times[-1])
    for i, agent in enumerate(agents):
        x_timeseries = path_timeseries[
            PATH_TO_AGENTS + (agent,)
            + path_to_x_variable][min_idx:max_idx]
        y_timeseries = path_timeseries[
            PATH_TO_AGENTS + (agent,)
            + path_to_y_variable][min_idx:max_idx]
        dead_timeseries = path_timeseries[
            PATH_TO_AGENTS + (agent,)
            + PATH_TO_DEAD][min_idx:max_idx]
        color_timeseries = tuple(
            DEAD_COLOR if dead else LIVE_COLOR
            for dead in dead_timeseries)
        ax.scatter(
            np.array(x_timeseries) * scaling,
            np.array(y_timeseries) * scaling,
            c=color_timeseries,
            label=agent,
            marker=MARKERS[i],
            alpha=ALPHA)
        ax.plot(
            np.array(x_timeseries) * scaling,
            np.array(y_timeseries) * scaling,
            color='black',
            linewidth=0.5)
    ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    for spine_name in ('top', 'right'):
        ax.spines[spine_name].set_visible(False)
    fig.tight_layout()
    return fig


def plot_expression_survival_dotplot(
    data, path_to_variable, xlabel, scaling=1, time_range=(0, 1)
):
    '''Create Expression Dotplot Colored by Survival

    Plot one dot for each cell along an axis to indicate that cell's
    average expression level for a specified protein. The dot color
    reflects whether the cell survived long enough to divide.

    Note that only the expression levels while the cell is alive are
    considered in the average.

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
            the time range to consider when calculating average
            expression level.

    Returns:
        plt.Figure: The finished figure.
    '''
    live_averages, dead_averages = calc_live_and_dead_averages(
        data, path_to_variable, time_range)
    fig, ax = plt.subplots(figsize=(6, 2))
    ax.scatter(
        np.array(list(live_averages.values())) * scaling,
        [0.1] * len(live_averages),
        label='Survive', color=LIVE_COLOR, alpha=ALPHA,
    )
    ax.scatter(
        np.array(list(dead_averages.values())) * scaling,
        [0.1] * len(dead_averages),
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


def calc_live_and_dead_averages(data, path_to_variable, time_range):
    expression_levels = dict()
    die = set()
    end_time = max(data.keys())
    for time, time_data in data.items():
        if (time < time_range[0] * end_time
                or time > time_range[1] * end_time):
            continue
        agents_data = get_in(time_data, PATH_TO_AGENTS)
        for agent, agent_data in agents_data.items():
            lst = expression_levels.setdefault(agent, [])
            value = get_in(agent_data, path_to_variable)
            if value is not None:
                lst.append(value)
            if get_in(agent_data, PATH_TO_DEAD, False):
                die.add(agent)
            # Only count values when cell is alive
            elif value is not None:
                lst.append(value)

    live_averages = {}
    dead_averages = {}
    for agent, levels in expression_levels.items():
        if not levels:
            continue
        if agent in die:
            dead_averages[agent] = np.mean(levels)
        else:
            live_averages[agent] = np.mean(levels)
    return live_averages, dead_averages
