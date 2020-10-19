'''Script for investigating why colony imports, then exports antibiotic

This import and export oscillation may not even occur anymore, which
this investigation will also reveal.
'''

import os
from typing import Iterable, List

from matplotlib import pyplot as plt
import numpy as np
from vivarium.core.experiment import get_in
from vivarium_cell.analysis.analyze import Analyzer
from vivarium_cell.plots.multibody_physics import plot_snapshots

from src.constants import FIELDS_PATH
from src.investigate_utils import (
    parse_args_retrieve_data_by_id,
    filter_raw_data_by_time,
)
from src.types import RawData


TIME_RANGE = (0.51, 1)
FIELDS = ['nitrocefin']
N_SNAPSHOTS = 35


def get_field_averages(raw_data: RawData, field: str) -> List[float]:
    '''Calculate the average concentration for a field over time.

    Args:
        raw_data: The raw simulation data.
        field: The field variable name.

    Returns:
        A list, ordered by time, of the average field concentration
        across all cells in the lattice at each timepoint.
    '''
    fields_timeseries = [
        get_in(raw_data[time], FIELDS_PATH + (field,))
        for time in sorted(raw_data.keys())
    ]
    field_averages = [
        np.mean(field)
        for field in fields_timeseries
    ]
    return field_averages


def plot_field_averages(
        raw_data: RawData, fields: Iterable[str]) -> plt.Figure:
    '''Plot the average concentration for fields over time.

    Args:
        raw_data: Raw simulation data.
        fields: The fields to plot.

    Returns:
        The figure.
    '''
    fig, ax = plt.subplots()
    times = sorted(raw_data.keys())
    for field in fields:
        averages = get_field_averages(raw_data, field)
        ax.plot(times, averages, label=field)
    return fig


def main() -> None:
    '''Run analyses'''
    data, environment_config, out_dir = parse_args_retrieve_data_by_id()
    data = filter_raw_data_by_time(data, TIME_RANGE)
    snapshots_data = Analyzer.format_data_for_snapshots(
        data, environment_config)
    plot_config = {
        'out_dir': out_dir,
        'filename': 'snapshots_during_pulse',
        'include_fields': FIELDS,
        'n_snapshots': N_SNAPSHOTS,
    }
    plot_snapshots(snapshots_data, plot_config)

    field_averages_fig = plot_field_averages(data, FIELDS)
    field_averages_fig.savefig(
        os.path.join(out_dir, 'field_avg_during_pulse.pdf'))


if __name__ == '__main__':
    main()
