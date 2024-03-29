'''Process Expression Data for Analysis'''

import argparse
import csv
import json
import os

from typing import Dict
import pandas as pd
from vivarium.library.topology import get_in

from src.db import add_connection_args, get_experiment_data
from src.types import Path, RawData
from src.constants import AGENTS_PATH, VOLUME_PATH


VOLUME_KEY = 'volume'


def raw_data_to_end_expression_table(
        raw_data: RawData,
        paths_dict: Dict[str, Path]) -> pd.DataFrame:
    '''Create a table of end expression levels from raw simulation data.

    Args:
        raw_data: Raw simulation data
        paths: Map from names to paths to protein counts. The names will
            be used as column headers in the returned table.

    Returns:
        Table with one column for each protein and one row for each
        agent. Each cell contains the protein concentration in that
        agent in the final simulation timepoint.
    '''
    end_data = raw_data[max(raw_data.keys())]
    expression_data: Dict = {
        name: []
        for name in paths_dict
    }
    expression_data[VOLUME_KEY] = []
    agents_data = get_in(end_data, AGENTS_PATH)
    for agent_data in agents_data.values():
        volume = get_in(agent_data, VOLUME_PATH, 0)
        expression_data[VOLUME_KEY].append(volume)
        for name, path in paths_dict.items():
            count = get_in(agent_data, path, 0)
            concentration = count / volume if volume else 0
            expression_data[name].append(concentration)
    return pd.DataFrame(expression_data)


def process_data(args: argparse.Namespace) -> None:
    '''Process expression data

    Args:
        args: Command-line arguments from argparse.
    '''
    if args.json_file:
        with open(args.json_path, 'r') as f:
            data = json.load(f)
    else:
        data, _ = get_experiment_data(args, args.experiment_id)
    with open(args.tagged_molecules, 'r') as f:
        reader = csv.reader(f)
        paths = [tuple(line) for line in reader]
    paths_dict = {
        path[-1]: path for path in paths
    }
    table = raw_data_to_end_expression_table(data, paths_dict)
    out_dir = os.path.dirname(args.output_csv)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)
    table.to_csv(args.output_csv)


def main() -> None:
    '''Main function that handles CLI arguments.'''
    parser = argparse.ArgumentParser()
    add_connection_args(parser)
    parser.add_argument(
        '-e', '--experiment_id',
        type=str,
        default=None,
        help='ID of experiment. Only needed to retrieve from database.'
    )
    parser.add_argument(
        '--json_file',
        type=str,
        default=None,
        help=(
            'JSON file with simulation output. If specified, no '
            'database is used.'
        ),
    )
    parser.add_argument(
        'tagged_molecules',
        type=str,
        default=None,
        help='Path to CSV file with tagged molecules paths.',
    )
    parser.add_argument(
        '-c', '--output_csv',
        type=str,
        default='expression.csv',
        help='Path to output CSV file.'
    )
    args = parser.parse_args()
    process_data(args)


if __name__ == '__main__':
    main()
