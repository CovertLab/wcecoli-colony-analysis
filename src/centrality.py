from typing import Sequence, Tuple

from matplotlib import pyplot as plt
import numpy as np
from vivarium.library.topology import get_in

from src.types import RawData
from src.investigate_utils import (
    split_raw_data_by_survival)
from src.constants import AGENTS_PATH, BOUNDS_PATH


LOCATION_PATH = ('boundary', 'location')

Location = Sequence[float]
Locations = Sequence[Location]


def get_survival_against_centrality_plot(
        data: RawData,
        fontsize: float = 20) -> Tuple[plt.Figure, dict]:
    '''Get a centrality box plot.

    Args:
        data: Simulation data.
        fontsize: Size of text on plot.

    Returns:
        Tuple of ``(fig, stats)`` with Matplotlib Figure object as
        ``fig`` and a dictionary of statistics as ``stats``.
    '''
    survive_locations, die_locations = extract_spatial_data(data)
    center = extract_center(data)

    fig, ax = plt.subplots()
    stats = plot_survival_against_centrality(
        survive_locations, die_locations, center, ax, fontsize)
    fig.tight_layout()
    return fig, stats


def extract_spatial_data(data: RawData) -> Tuple[Locations, Locations]:
    '''Get agent locations at the end of a simulation split by survival.

    Args:
        data: Simulation data.

    Returns:
        Tuple ``(survive_locations, die_locations)``.
    '''
    end_time = max(data.keys())
    data_filtered = RawData({
        end_time: data[end_time]
    })
    survive_data, die_data = split_raw_data_by_survival(data_filtered)
    survive_locations = [
        get_in(agent_data, LOCATION_PATH)
        for agent_data in get_in(
            survive_data[end_time],  # pylint: disable=unsubscriptable-object
            AGENTS_PATH).values()
    ]
    die_locations = [
        get_in(agent_data, LOCATION_PATH)
        for agent_data in get_in(
            die_data[end_time], AGENTS_PATH).values()
    ]
    return survive_locations, die_locations


def extract_center(data: RawData) -> Location:
    '''Find the environment center from simulation data.

    Args:
        data: Simulation raw data (top-level keys are times).

    Returns:
        Coordinates of the environment center.
    '''
    end_time = max(data.keys())
    bounds = get_in(data[end_time], BOUNDS_PATH)
    x, y = bounds
    return x / 2, y / 2


def plot_survival_against_centrality(
        survive_locations: Locations, die_locations: Locations,
        center: Location, ax: plt.Axes,
        fontsize: float = 20) -> dict:
    '''Create box plot of agent distances from center split by survival.

    Args:
        survive_locations: Locations of agents that survived to divide
            or until the end of the simulation.
        die_locations: Locations of agents that died.
        center: Coordinates of the environment center.
        ax: Axes to plot on.
        fontsize: Size of font for plot text.

    Returns:
        Statistics dictionary with all the distances of surviving agents
        to the center (under key ``survive_distances``) and all the
        distances of dying agents (under key ``die_distances``).
    '''
    survive_array = np.array(survive_locations)
    die_array = np.array(die_locations)
    center_array = np.array(center)

    to_plot = []
    if len(survive_array) != 0:
        to_plot.append(np.linalg.norm(  # type: ignore
            survive_array - center_array, ord=2, axis=1)
        )
    else:
        to_plot.append([])
    if len(die_array) != 0:
        to_plot.append(np.linalg.norm(  # type: ignore
            die_array - center_array, ord=2, axis=1)
        )
    else:
        to_plot.append([])

    median_props = {'color': 'black'}
    boxplot_labels = ('Survive', 'Die')
    ax.boxplot(  # type: ignore
        to_plot, labels=boxplot_labels, medianprops=median_props,
        showfliers=False)
    ax.tick_params(  # type: ignore
        axis='both', which='major', labelsize=fontsize)
    for i, y_values in enumerate(to_plot):
        ax.scatter(
            np.random.normal(i + 1, 0.04, size=len(y_values)),
            y_values,
            c='black' if boxplot_labels[i] == 'Die' else 'green',
            alpha=0.2,
        )
    ax.set_ylabel(  # type: ignore
        'Distance from Center ($\\mu m$)', fontsize=fontsize)
    for spine_name in ('top', 'right'):
        ax.spines[spine_name].set_visible(False)  # type: ignore

    stats = {
        'survive_distances': to_plot[0],
        'die_distances': to_plot[1],
    }
    return stats
