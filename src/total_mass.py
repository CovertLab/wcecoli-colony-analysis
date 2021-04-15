'''Colony Total Mass'''

from typing import Dict, List, Tuple, Iterable, cast

from matplotlib import pyplot as plt
import numpy as np
from vivarium.core.experiment import get_in

from src.constants import AGENTS_PATH, MASS_PATH
from src.types import RawData


def get_total_mass(agents_data: Dict) -> float:
    '''Get the total mass of a set of agents.

    Args:
        agents_data: The simulation data of the store containing all the
            agents.

    Returns:
        The total mass of the agents.
    '''
    total_mass = 0.
    for agent_data in agents_data.values():
        mass = get_in(agent_data, MASS_PATH)
        total_mass += mass
    return total_mass


def get_total_mass_timeseries(data: RawData) -> List[float]:
    '''Get a timeseries of the total mass of a simulation.

    Args:
        data: Data from the simulation.

    Returns:
        A list of the total cell mass in the simulation over time.
    '''
    times = sorted(data.keys())
    mass_timeseries = []
    for time in times:
        agents_data = get_in(data[time], AGENTS_PATH)
        mass_timeseries.append(get_total_mass(agents_data))
    return mass_timeseries


def get_total_mass_plot(
        datasets: Dict[str, List[RawData]],
        colors: List[str],
        fontsize: float = 36,
        vlines: Iterable[Tuple[float, float, str, str]] = tuple(),
        ) -> Tuple[plt.Figure, dict]:
    '''Plot the total masses of colonies from groups of simulations.

    Each group's total mass over time is plotted as a curve on the
    resulting figure.

    Args:
        datasets: Map from the label to associate with a group of
            simulations to a list of the datasets in that group.
        colors: Map from a group label to the color to show that group's
            data in.
        fontsize: Size of all text on figure.
        vlines: Tuple of vertical line specifiers. Each specifier is a
            tuple of the line position, label position as fraction of x
            range, color, and label.

    Returns:
        A tuple of the figure and a dictionary that maps from group
        label to tuples of the first, second, and third quartiles of
        that group's data.
    '''
    fig, ax = plt.subplots()
    quartiles = {}
    for i, (label, replicates) in enumerate(datasets.items()):
        filtered_replicates = []
        for replicate in replicates:
            # Exclude first timepoint, which is often wrong
            filtered = RawData({
                key: val for key, val in replicate.items()
                if key != min(replicate.keys())
            })
            filtered_replicates.append(filtered)
        label_quartiles = plot_total_mass(
            filtered_replicates, ax, label, colors[i], fontsize)
        quartiles[label] = label_quartiles
    for x, label_x, vline_color, vline_label in vlines:
        ax.axvline(  # type: ignore
            x / 60 / 60, color=vline_color, linestyle='--')
        ax.text(  # type: ignore
            label_x, 0.95, vline_label, fontsize=fontsize,
            transform=ax.transAxes)  # type: ignore
    ax.set_ylabel(  # type: ignore
        'Total Cell Mass (fg)', fontsize=fontsize)
    ax.set_xlabel('Time (hr)', fontsize=fontsize)  # type: ignore
    for spine_name in ('top', 'right'):
        ax.spines[spine_name].set_visible(False)  # type: ignore
    fig.tight_layout()
    return fig, quartiles


def plot_total_mass(
        replicates: List[RawData],
        ax: plt.Axes,
        label: str = '',
        color: str = 'black',
        fontsize: float = 36,
        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''Plot total mass data on an existing set of axes.

    Plots the median surrounded by a translucent band indicating the
    IQR. The median and IQR are computed over the replicates.

    Args:
        replicates: A list of the raw data dictionary for each
            replicate.
        ax: The axes to plot on.
        label: Label to associate with the curve showing the median
            total mass.
        color: Color of median curve and IQR band.
        fontsize: Size of all text in figure.

    Returns:
        A tuple of a numpy array for each quartile.
    '''
    times = sorted(replicates[0].keys())
    mass_timeseries = []
    for replicate in replicates:
        assert sorted(replicate.keys()) == times
        replicate_timeseries = get_total_mass_timeseries(replicate)
        mass_timeseries.append(replicate_timeseries)
    mass_matrix = np.array(mass_timeseries)
    median = np.median(mass_matrix, axis=0)
    q25, q75 = np.percentile(mass_matrix, [25, 75], axis=0)
    times_hours = tuple(time / 60 / 60 for time in times)
    if label:
        ax.semilogy(  # type: ignore
            times_hours, median, label=label, color=color)
        ax.fill_between(  # type: ignore
            times_hours, q25, q75, color=color, alpha=0.2,
            edgecolor='none')
        ax.legend(  # type: ignore
            prop={'size': fontsize}, frameon=False)
    else:
        ax.semilogy(  # type: ignore
            times_hours, mass_timeseries, color=color)
        ax.fill_between(  # type: ignore
            times_hours, q25, q75, color=color, alpha=0.2,
            edgecolor='none')
    for tick_type in ('major', 'minor'):
        ax.tick_params(  # type: ignore
            axis='both', which=tick_type, labelsize=fontsize)
    return cast(np.ndarray, q25), median, cast(np.ndarray, q75)
