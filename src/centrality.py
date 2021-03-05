from typing import Sequence, Tuple

from matplotlib import pyplot as plt
import numpy as np
from vivarium.core.experiment import get_in

from src.types import RawData
from src.investigate_utils import (
    split_raw_data_by_survival, filter_raw_data_by_time)
from src.constants import AGENTS_PATH, BOUNDS_PATH


LOCATION_PATH = ('boundary', 'location')

Location = Sequence[float]
Locations = Sequence[Location]


def get_survival_against_centrality_plot(data: RawData) -> plt.Figure:
    survive_locations, die_locations = extract_spatial_data(data)
    center = extract_center(data)

    fig, ax = plt.subplots()
    plot_survival_against_centrality(
        survive_locations, die_locations, center, ax)
    return fig


def extract_spatial_data(data: RawData) -> Tuple[Locations, Locations]:
    end_time = max(data.keys())
    data_filtered = RawData({
        end_time: data[end_time]
    })
    survive_data, die_data = split_raw_data_by_survival(data_filtered)
    survive_locations = [
        get_in(agent_data, LOCATION_PATH)
        for agent_data in get_in(
            survive_data[end_time], AGENTS_PATH).values()
    ]
    die_locations = [
        get_in(agent_data, LOCATION_PATH)
        for agent_data in get_in(
            die_data[end_time], AGENTS_PATH).values()
    ]
    return survive_locations, die_locations


def extract_center(data: RawData) -> Location:
    bounds = get_in(data, BOUNDS_PATH)
    x, y = bounds
    return x / 2, y / 2


def plot_survival_against_centrality(
        survive_locations: Locations, die_locations: Locations,
        center: Location, ax: plt.Axes) -> None:
    survive_array = np.array(survive_locations)  # type: ignore
    die_array = np.array(die_locations)  # type: ignore
    center_array = np.array(center)  # type: ignore

    survive_distances = np.linalg.norm(  # type: ignore
        survive_array - center_array, ord=2, axis=1)
    die_distances = np.linalg.norm(  # type: ignore
        die_array - center_array, ord=2, axis=1)

    ax.boxplot(  # type: ignore
        [survive_distances, die_distances], labels=['Survive', 'Die'])
    for i, y_values in enumerate((survive_distances, die_distances)):
        ax.scatter(
            np.random.normal(i + 1, 0.04, size=len(y_values)),
            y_values,
            c='black', alpha=0.2,
        )
    ax.set_ylabel('Euclidian Distance from Environment Center')
