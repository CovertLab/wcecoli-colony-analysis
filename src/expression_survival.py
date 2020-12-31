'''
======================================
Expression Dotplot Colored by Survival
======================================

Adapted from vivarium-cell.
'''
from matplotlib import pyplot as plt
import numpy as np

from vivarium.core.experiment import get_in


PATH_TO_AGENTS = ('agents',)
PATH_TO_DEAD = ('boundary', 'dead')
LIVE_COLOR = 'red'
DEAD_COLOR = 'black'
ALPHA = 0.5


def plot_expression_survival(
    data, path_to_x_variable, path_to_y_variable, xlabel, ylabel, time_range=(0, 1)
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
        time_range (tuple): Tuple of two :py:class:`float`s that are
            fractions of the total simulated time period. These
            fractions indicate the start and end points (inclusive) of
            the time range to consider when calculating average
            expression level.

    Returns:
        plt.Figure: The finished figure.
    '''
    live_averages_x, dead_averages_x = calc_live_and_dead_averages(
        data, path_to_x_variable, time_range)
    live_averages_y, dead_averages_y = calc_live_and_dead_averages(
        data, path_to_y_variable, time_range)
    fig, ax = plt.subplots(figsize=(6, 2))
    ax.scatter(
        live_averages_x, live_averages_y,
        label='Survive', color=LIVE_COLOR, alpha=ALPHA,
    )
    ax.scatter(
        dead_averages_y, dead_averages_y,
        label='Die', color=DEAD_COLOR, alpha=ALPHA,
    )
    ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    for spine_name in ('top', 'right'):
        ax.spines[spine_name].set_visible(False)
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

    live_averages = []
    dead_averages = []
    for agent, levels in expression_levels.items():
        if not levels:
            continue
        if agent in die:
            dead_averages.append(np.mean(levels))
        else:
            live_averages.append(np.mean(levels))
    return live_averages, dead_averages