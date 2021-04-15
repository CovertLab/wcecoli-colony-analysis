'''
======================================
Expression Dotplot Colored by Survival
======================================

Adapted from vivarium-cell.
'''
from typing import Sequence, Iterable, Tuple, List, Dict

from matplotlib import pyplot as plt
import numpy as np
from vivarium.core.experiment import get_in
from vivarium.core.emitter import path_timeseries_from_data

from src.investigate_utils import filter_raw_data_by_time
from src.types import RawData, Path


PATH_TO_AGENTS = ('agents',)
PATH_TO_DEAD = ('boundary', 'dead')
LIVE_COLOR = 'green'
DEAD_COLOR = 'black'
ALPHA = 0.5


def _get_final_live_agents(
        data: RawData,
        time_range: Tuple[float, float] = (0, 1),
        ) -> List[str]:
    data = filter_raw_data_by_time(data, time_range)
    max_time = max(data.keys())
    agents = []
    agents_data = get_in(
        # Pylint doesn't recognize that the RawData NewType is a dict
        data[max_time], # pylint: disable=unsubscriptable-object
        PATH_TO_AGENTS,
    )
    for agent, agent_data in agents_data.items():
        dead = get_in(agent_data, PATH_TO_DEAD)
        if not dead:
            agents.append(agent)
    return agents


def plot_expression_survival(
        data: RawData,
        path_to_x_variable: Path,
        path_to_y_variable: Path,
        xlabel: str,
        ylabel: str,
        boundary_x: Sequence[float],
        boundary_y: Sequence[float],
        boundary_error: Sequence[float],
        boundary_color: str = 'black',
        scaling: float = 1,
        time_range: Tuple[float, float] = (0, 1),
        label_agents: bool = False,
        plot_agents: Iterable[str] = tuple(),
        fontsize: float = 36,
        dead_trace_agents: Iterable[str] = tuple(),
        agents_for_phylogeny_trace: Iterable[str] = tuple(),
        ) -> plt.Figure:
    '''Create Expression Scatterplot Colored by Survival

    Plot one dot for each cell along an axis to indicate that cell's
    final expression level for a specified protein. The dot color
    reflects whether the cell survived long enough to divide.

    Parameters:
        data: The raw data emitted from the simulation.
        path_to_x_variable: Path from the agent root to the
            variable that holds the protein's expression level. We do
            not adjust for cell volume, so this should be a
            concentration. This protein will be plotted on the x axis.
        path_to_y_variable: Path from the agent root to the
            variable that holds the protein's expression level. We do
            not adjust for cell volume, so this should be a
            concentration. This protein will be plotted on the y axis.
        xlabel: Label for x-axis.
        xlabel: Label for y-axis.
        boundary_x: X-values of the boundary identified
            numerically from the antibiotic model.
        boundary_y: Y-values of the boundary identified
            numerically from the antibiotic model.
        boundary_error: Precision of the Y-value
            predictions. This is the distance along the y axis between
            the known dead point and known live point closest to the
            prediction. The predicted y-value will be in the middle of
            this range, so an error band can be drawn of width
            boundary_error centered at boundary_y.
        boundary_color: Color of boundary.
        scaling: Coefficient to multiply all data by. This is
            intended to be used for changing the units plotted.
        time_range: Tuple of two :py:class:`float`s that are
            fractions of the total simulated time period. These
            fractions indicate the start and end points (inclusive) of
            the time range to consider.
        label_agents: Whether to label each point with the agent
            ID.
        plot_agents: The agent IDs of the agents to plot. By
            default, all agents are plotted.
        fontsize: Text size for entire figure.
        dead_trace_agents: The agent IDs of the agents to
            plot traces for. By default, no traces are shown.
        agents_for_phylogeny_trace: Agent IDs for the agents
            whose phylogenies will be traced.

    Returns:
        The finished figure.
    '''
    live_finals_x, dead_finals_x = _calc_live_and_dead_finals(
        data, path_to_x_variable, time_range, plot_agents)
    live_finals_y, dead_finals_y = _calc_live_and_dead_finals(
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

    if live_finals_x:
        ax.scatter(  # type: ignore
            np.array(list(live_finals_x.values())) * scaling,
            np.array(list(live_finals_y.values())) * scaling,
            label='Concentrations at Division (Cell Survives)',
            color=LIVE_COLOR, alpha=ALPHA,
        )
    if dead_finals_x:
        ax.scatter(  # type: ignore
            np.array(list(dead_finals_x.values())) * scaling,
            np.array(list(dead_finals_y.values())) * scaling,
            label='Concentrations at Death', color=DEAD_COLOR,
            alpha=ALPHA,
        )
    if label_agents:
        for agent in live_finals_x:
            x = live_finals_x[agent] * scaling
            y = live_finals_y[agent] * scaling
            ax.annotate(agent, (x, y), size=0.1)  # type: ignore
        for agent in dead_finals_x:
            x = dead_finals_x[agent] * scaling
            y = dead_finals_y[agent] * scaling
            ax.annotate(agent, (x, y), size=0.1)  # type: ignore
    plot_expression_survival_death_traces(
        ax, data, path_to_x_variable, path_to_y_variable, scaling,
        time_range, dead_trace_agents, DEAD_COLOR)
    plot_expression_survival_lineage_traces(
        ax, data, path_to_x_variable, path_to_y_variable, scaling,
        time_range, agents_for_phylogeny_trace, LIVE_COLOR, ALPHA)
    finals = list(live_finals_x.values()) + list(
        dead_finals_x.values())
    plot_expression_survival_boundary(
        ax, boundary_x, boundary_y, boundary_error, finals, scaling,
        boundary_color)
    ax.legend(  # type: ignore
        bbox_to_anchor=(0.5, 1.05), loc='lower center',
        frameon=False, prop={'size': fontsize})
    ax.set_xlabel(xlabel, fontsize=fontsize)  # type: ignore
    ax.set_ylabel(ylabel, fontsize=fontsize)  # type: ignore
    ax.tick_params(  # type: ignore
        axis='both', which='major', labelsize=fontsize)
    for spine_name in ('top', 'right'):
        ax.spines[spine_name].set_visible(False)  # type: ignore
    fig.tight_layout()
    return fig


def plot_expression_survival_boundary(
        ax: plt.Axes,
        boundary_x: Sequence[float],
        boundary_y: Sequence[float],
        boundary_error: Sequence[float],
        finals: Sequence[float],
        scaling: float = 1,
        boundary_color: str = 'black'
        ) -> None:
    '''Plot the life-death boundary on an expression-survival plot.

    Args:
        ax: Axes to plot on.
        boundary_x: X-values of the boundary identified
            numerically from the antibiotic model.
        boundary_y: Y-values of the boundary identified
            numerically from the antibiotic model.
        boundary_error: Precision of the Y-value
            predictions. This is the distance along the y axis between
            the known dead point and known live point closest to the
            prediction. The predicted y-value will be in the middle of
            this range, so an error band can be drawn of width
            boundary_error centered at boundary_y.
        finals: All final cell concentrations plotted on
            the figure. This is used to ensure the boundary line spans
            the figure.
        scaling: Coefficient to multiply all data by. This is
            intended to be used for changing the units plotted.
        boundary_color: Color of boundary.
    '''
    boundary_x_arr = np.array(boundary_x)
    boundary_y_arr = np.array(boundary_y)
    boundary_error_arr = np.array(boundary_error)
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
    ax.plot(  # type: ignore
        boundary_x_arr * scaling, boundary_y_arr * scaling,
        c=boundary_color,
        label='Boundary Predicted by Antibiotics Model')
    ax.fill_between(  # type: ignore
        boundary_x_arr * scaling,
        boundary_y_arr * scaling - boundary_error_arr * scaling / 2,
        boundary_y_arr * scaling + boundary_error_arr * scaling / 2,
        color=boundary_color, alpha=0.2)



def plot_expression_survival_death_traces(
        ax: plt.Axes,
        data: RawData,
        path_to_x_variable: Path,
        path_to_y_variable: Path,
        scaling: float = 1,
        time_range: Tuple[float, float] = (0, 1),
        dead_agents: Iterable[str] = tuple(),
        dead_trace_color: str = 'black',
        ) -> None:
    '''Create Expression Traces for Dead Cells

    Parameters:
        ax: Axes to plot on.
        data: The raw data emitted from the simulation.
        path_to_x_variable: Path from the agent root to the
            variable that holds the protein's expression level. We do
            not adjust for cell volume, so this should be a
            concentration. This protein will be plotted on the x axis.
        path_to_y_variable: Path from the agent root to the
            variable that holds the protein's expression level. We do
            not adjust for cell volume, so this should be a
            concentration. This protein will be plotted on the y axis.
        scaling: Coefficient to multiply all data by. This is
            intended to be used for changing the units plotted.
        time_range: Tuple of two :py:class:`float`s that are
            fractions of the total simulated time period. These
            fractions indicate the start and end points (inclusive) of
            the time range to consider.
        dead_agents: The agent IDs of the agents to plot
            traces for. These agents should die.
        dead_trace_color: Color of trace line for dead cells.
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
        ax.plot(  # type: ignore
            np.array(x_timeseries) * scaling,
            np.array(y_timeseries) * scaling,
            color=dead_trace_color,
            linewidth=1,
            label='Agent path until death' if i == 0 else '',
        )


def plot_expression_survival_lineage_traces(
        ax: plt.Axes,
        data: RawData,
        path_to_x_variable: Path,
        path_to_y_variable: Path,
        scaling: float = 1,
        time_range: Tuple[float, float] = (0, 1),
        agents_for_phylogeny_trace: Iterable[str] = tuple(),
        phylogeny_trace_color: str = 'green',
        alpha: float = 1,
        ) -> None:
    '''Create expression traces for a lineage of cells.

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
        agents_for_phylogeny_trace (Iterable): Agent IDs for the agents
            whose phylogenies will be traced.
        phylogeny_trace_color (str): Color of trace line for phylogeny.
        alpha (float): Transparency for starting point.
    '''
    data = filter_raw_data_by_time(data, time_range)
    path_timeseries = path_timeseries_from_data(data)

    # Plot phylogeny traces
    plotted_solid = False
    plotted_dashed = False
    for agent in agents_for_phylogeny_trace:
        last_end_point: Tuple[float, ...] = tuple()
        for i in range(len(agent) + 1):
            ancestor = agent[:i]
            x_path = PATH_TO_AGENTS + (ancestor,) + path_to_x_variable
            y_path = PATH_TO_AGENTS + (ancestor,) + path_to_y_variable
            if x_path in path_timeseries:
                # else ancestor does not exist
                x_timeseries = path_timeseries[x_path]
                y_timeseries = path_timeseries[y_path]
                ax.plot(  # type: ignore
                    np.array(x_timeseries) * scaling,
                    np.array(y_timeseries) * scaling,
                    color=phylogeny_trace_color,
                    linewidth=1,
                    label=(
                        'Agent path until division'
                        if not plotted_solid else ''
                    ),
                )
                if not plotted_solid:
                    # This is the first agent in the lineage
                    ax.scatter(  # type: ignore
                        x_timeseries[0] * scaling,
                        y_timeseries[0] * scaling,
                        color=phylogeny_trace_color,
                        marker='s',
                        label='Lineage start',
                        alpha=alpha,
                    )
                plotted_solid = True
                if last_end_point:
                    ax.plot(  # type: ignore
                        [
                            last_end_point[0] * scaling,
                            x_timeseries[0] * scaling
                        ],
                        [
                            last_end_point[1] * scaling,
                            y_timeseries[0] * scaling
                        ],
                        color=phylogeny_trace_color,
                        linewidth=1,
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
        data: RawData,
        path_to_variable: Path,
        xlabel: str,
        scaling: float = 1,
        time_range: Tuple[float, float] = (0, 1),
        fontsize: float = 36,
    ) -> Tuple[plt.Figure, dict]:
    '''Create Expression Dotplot Colored by Survival

    Plot one dot for each cell along an axis to indicate that cell's
    final expression level for a specified protein. The dot color
    reflects whether the cell survived long enough to divide.

    Note that only the expression levels while the cell is alive are
    considered.

    Parameters:
        data: The raw data emitted from the simulation.
        path_to_variable: Path from the agent root to the
            variable that holds the protein's expression level. We do
            not adjust for cell volume, so this should be a
            concentration.
        xlabel: Label for x-axis.
        scaling: Coefficient to multiply all data by. This is
            intended to be used for changing the units plotted.
        time_range: Tuple of two :py:class:`float`s that are
            fractions of the total simulated time period. These
            fractions indicate the start and end points (inclusive) of
            the time range to consider.
        fontsize: size of all text in figure.

    Returns:
        Tuple of The finished figure and a statistics dictionary with
        the final concentrations of the live cells (under key ``live``)
        and the dead cells (under key ``dead``).
    '''
    live_finals, dead_finals = _calc_live_and_dead_finals(
        data, path_to_variable, time_range)
    fig, ax = plt.subplots(figsize=(6, 2))
    ax.scatter(  # type: ignore
        np.array(list(live_finals.values())) * scaling,
        [0.1] * len(live_finals),
        label='Survive', color=LIVE_COLOR, alpha=ALPHA,
    )
    ax.scatter(  # type: ignore
        np.array(list(dead_finals.values())) * scaling,
        [0.1] * len(dead_finals),
        label='Die', color=DEAD_COLOR, alpha=ALPHA,
    )
    ax.legend(prop={'size': fontsize}, frameon=False)  # type: ignore
    ax.set_xlabel(xlabel, fontsize=fontsize)  # type: ignore
    ax.set_ylim([0, 1.25])  # type: ignore
    ax.get_yaxis().set_visible(False)  # type: ignore
    for spine_name in ('left', 'top', 'right'):
        ax.spines[spine_name].set_visible(False)  # type: ignore
    ax.spines['bottom'].set_position('zero')  # type: ignore
    fig.tight_layout()
    stats = {
        'live': list(live_finals.values()),
        'dead': list(dead_finals.values()),
    }
    return fig, stats


def _calc_live_and_dead_finals(
        data: RawData,
        path_to_variable: Path,
        time_range: Tuple[float, float],
        agents: Iterable[str] = tuple(),
        ) -> Tuple[Dict[str, float], Dict[str, float]]:
    values: Dict[str, Dict[float, float]] = {}
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
