from typing import Dict, List

from matplotlib import pyplot as plt
import numpy as np
from vivarium.core.experiment import get_in

from src.constants import AGENTS_PATH, MASS_PATH
from src.types import RawData


def get_total_mass(agents_data: Dict) -> float:
    total_mass = 0.
    for agent_data in agents_data.values():
        mass = get_in(agent_data, MASS_PATH)
        total_mass += mass
    return total_mass


def get_total_mass_timeseries(data: RawData) -> List[float]:
    times = sorted(data.keys())
    mass_timeseries = []
    for time in times:
        agents_data = get_in(data[time], AGENTS_PATH)
        mass_timeseries.append(get_total_mass(agents_data))
    return mass_timeseries


def get_total_mass_plot(
        datasets: Dict[str, List[RawData]],
        colors: List[str]) -> plt.Figure:
    fig, ax = plt.subplots()
    for i, (label, replicates) in enumerate(datasets.items()):
        filtered_replicates = []
        for replicate in replicates:
            # Exclude first timepoint, which is often wrong
            filtered = RawData({
                key: val for key, val in replicate.items()
                if key != min(replicate.keys())
            })
            filtered_replicates.append(filtered)
        plot_total_mass(filtered_replicates, ax, label, colors[i])
    ax.set_ylabel('Total Cell Mass (fg)')
    ax.set_xlabel('Time (s)')
    fig.tight_layout()
    return fig


def plot_total_mass(
        replicates: List[RawData], ax: plt.Axes, label='',
        color='black') -> None:
    times = sorted(replicates[0].keys())
    mass_timeseries = []
    for replicate in replicates:
        assert sorted(replicate.keys()) == times
        replicate_timeseries = get_total_mass_timeseries(replicate)
        mass_timeseries.append(replicate_timeseries)
        ax.semilogy(
            times, replicate_timeseries, color=color, linestyle='--')
    mass_matrix = np.array(mass_timeseries)
    median = np.median(mass_matrix, axis=0)
    q25, q75 = np.percentile(mass_matrix, [25, 75], axis=0)
    if label:
        ax.semilogy(times, median, label=label, color=color)
        ax.fill_between(times, q25, q75, color=color, alpha=0.2)
        ax.legend()
    else:
        ax.semilogy(times, mass_timeseries, color=color)
        ax.fill_between(times, q25, q75, color=color, alpha=0.2)
