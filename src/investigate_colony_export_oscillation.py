'''Script for investigating why colony imports, then exports antibiotic

This import and export oscillation may not even occur anymore, which
this investigation will also reveal.
'''

from vivarium_cell.analysis.analyze import Analyzer
from vivarium_cell.plots.multibody_physics import plot_snapshots

from src.investigate_utils import (
    parse_args_retrieve_data_by_id,
    filter_raw_data_by_time,
)


TIME_RANGE = (0.51, 1)
FIELDS = ['nitrocefin']
N_SNAPSHOTS = 50


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


if __name__ == '__main__':
    main()
