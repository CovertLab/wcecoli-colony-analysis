'''Generate video of simulation.

For usage information, run:
    python make_video.py -h
'''
import argparse
import os
from typing import List

from matplotlib import rcParams  # type: ignore
from vivarium_cell.analysis.analyze import Analyzer

from src.constants import OUT_DIR
from src.snapshots_video import make_tags_video, make_snapshots_video
from src.make_figures import get_experiment_data, TAG_PATH_NAME_MAP


FIG_OUT_DIR = os.path.join(OUT_DIR, 'figs')
VIDEO_CONFIGS = (
    {
        'name': 'glucose_consumption_basal',
        'experiment_id': '20201119.150828',
        'field': 'GLC',
        'video_type': 'snapshots',
        'field_added': '0',
    },
    {
        'name': 'glucose_consumption_anaerobic',
        'experiment_id': '20201221.194828',
        'field': 'GLC',
        'video_type': 'snapshots',
        'field_added': '0',
    },
    {
        'name': 'antibiotics',
        'experiment_id': '20210329.155953',
        'field': 'nitrocefin',
        'video_type': 'snapshots',
        'field_added': '0.5',
    },
    {
        'name': 'AmpC',
        'experiment_id': '20201119.150828',
        'tag': (
            'boundary', 'bulk_molecules_report', 'EG10040-MONOMER[p]'),
        'video_type': 'tags',
    },
    {
        'name': 'AcrABTolC',
        'experiment_id': '20201119.150828',
        'tag': (
            'boundary', 'bulk_molecules_report', 'TRANS-CPLX-201[s]'),
        'video_type': 'tags',
    },
)


def main() -> None:
    '''Generate all figures.'''
    rcParams['font.sans-serif'] = ['Arial']
    rcParams['font.family'] = ['sans-serif']
    if not os.path.exists(FIG_OUT_DIR):
        os.makedirs(FIG_OUT_DIR)
    parser = argparse.ArgumentParser(
        description=(
            'Generate selected figures and associated stats from '
            'simulation data.'
        ),
    )
    Analyzer.add_connection_args(parser)
    parser.add_argument(
        '--data_path',
        default='',
        help='Folder of JSON files to read data from instead of Mongo',
    )
    args = parser.parse_args()
    for config  in VIDEO_CONFIGS:
        data, environment_config = get_experiment_data(
            args, config['experiment_id'])
        if config['video_type'] == 'snapshots':
            make_snapshots_video(
                data, environment_config, [config['field']],
                filename=config['name'], out_dir=FIG_OUT_DIR,
                field_added=float(config['field_added']), xlim=(10, 40),
                ylim=(10, 40),
            )
        elif config['video_type'] == 'tags':
            make_tags_video(
                data, environment_config, [config['tag']],
                TAG_PATH_NAME_MAP, filename=config['name'],
                out_dir=FIG_OUT_DIR, xlim=(10, 40), ylim=(10, 40),
            )
        else:
            raise ValueError(
                'Unknown video type {}'.format(config['video_type']))


if __name__ == '__main__':
    main()
