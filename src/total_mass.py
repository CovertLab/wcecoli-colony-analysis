from typing import Dict, List

from matplotlib import pyplot as plt
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


def get_total_mass_plot(datasets: Dict[str, RawData]) -> plt.Figure:
    fig, ax = plt.subplots()
    for label, dataset in datasets.items():
        # Exclude first timepoint, which is often wrong
        filtered = RawData({
            key: val for key, val in dataset.items()
            if key != min(dataset.keys())
        })
        plot_total_mass(filtered, ax, label)
    ax.set_ylabel('Total Cell Mass (fg)')
    ax.set_xlabel('Time (s)')
    fig.tight_layout()
    return fig


def plot_total_mass(data: RawData, ax: plt.Axes, label='') -> None:
    times = sorted(data.keys())
    mass_timeseries = get_total_mass_timeseries(data)
    if label:
        ax.semilogy(times, mass_timeseries, label=label)
        ax.legend()
    else:
        ax.semilogy(times, mass_timeseries)
