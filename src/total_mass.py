'''Colony Total Mass'''

from typing import Dict, List, Tuple, cast

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
        colors: List[str]) -> Tuple[plt.Figure, dict]:
    '''Plot the total masses of colonies from groups of simulations.

    Each group's total mass over time is plotted as a curve on the
    resulting figure.

    Args:
        datasets: Map from the label to associate with a group of
            simulations to a list of the datasets in that group.
        colors: Map from a group label to the color to show that group's
            data in.

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
            filtered_replicates, ax, label, colors[i])
        quartiles[label] = label_quartiles
    ax.set_ylabel('Total Cell Mass (fg)')
    ax.set_xlabel('Time (s)')
    fig.tight_layout()
    return fig, quartiles


def plot_total_mass(
        replicates: List[RawData], ax: plt.Axes, label='',
        color='black') -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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
    if label:
        ax.semilogy(  # type: ignore
            times, median, label=label, color=color)
        ax.fill_between(  # type: ignore
            times, q25, q75, color=color, alpha=0.2, edgecolor='none')
        ax.legend()
    else:
        ax.semilogy(times, mass_timeseries, color=color)  # type: ignore
        ax.fill_between(  # type: ignore
            times, q25, q75, color=color, alpha=0.2, edgecolor='none')
    return cast(np.ndarray, q25), median, cast(np.ndarray, q75)
