import argparse
import json
import os

from vivarium.core.emitter import (
    get_local_client,
    data_from_database,
)

from src.types import RawData, EnvironmentConfig, DataTuple


def get_experiment_data(
        args: argparse.Namespace,
        experiment_id: str,
        ) -> DataTuple:
    '''Get simulation data for an experiment.

    If ``args.data_path`` is set, retrieve the experiment data from a
    JSON file named ``<experiment_id>.json`` under ``args.data_path``.
    Otherwise, retrieve the data from MongoDB.

    Args:
        args: Parsed CLI args.
        experiment_id: ID of experiment.

    Returns: Tuple of simulation data and environment config.
    '''
    if args.data_path:
        path = os.path.join(
            args.data_path, '{}.json'.format(experiment_id))
        with open(path, 'r') as f:
            loaded_file = json.load(f)
            data = RawData({
                float(time): value
                for time, value in loaded_file['data'].items()
            })
            config = EnvironmentConfig(
                loaded_file['environment_config'])
            return data, config
    client = get_local_client(
        args.host, args.port, args.database_name)
    return data_from_database(experiment_id, client)


def add_connection_args(parser: argparse.ArgumentParser) -> None:
    '''Add port, host, and db name args to an argument parser.'''
    parser.add_argument(
        '--port', '-p',
        default=27017,
        type=int,
        help=(
            'Port at which to access local mongoDB instance. '
            'Defaults to "27017".'
        ),
    )
    parser.add_argument(
        '--host', '-o',
        default='localhost',
        type=str,
        help=(
            'Host at which to access local mongoDB instance. '
            'Defaults to "localhost".'
        ),
    )
    parser.add_argument(
        '--database_name', '-d',
        default='simulations',
        type=str,
        help=(
            'Name of database on local mongoDB instance to read from. '
            'Defaults to "simulations".'
        )
    )


def format_data_for_snapshots(
        data: RawData,
        environment_config: EnvironmentConfig
        ) -> dict:
    '''Transform data into a form suitabe for generating snapshots.

    Args:
        data: The raw simulation data to transform.
        environment_config: Environment configuration.

    Returns:
        A dictionary with keys ``agents``, ``fields``, and ``config``
        suitable for passing to a snapshot plotting function.
    '''
    snapshots_data = {
        'agents': {
            time: timepoint['agents']
            for time, timepoint in data.items()
        },
        'fields': {
            time: timepoint['fields']
            for time, timepoint in data.items()
        },
        'config': environment_config,
    }
    return snapshots_data


def format_data_for_tags(
        data: RawData,
        environment_config: EnvironmentConfig
        ) -> dict:
    '''Transform data into a form suitabe for generating tag plots.

    Args:
        data: The raw simulation data to transform.
        environment_config: Environment configuration.

    Returns:
        A dictionary with keys ``agents`` and ``config`` suitable for
        passing to a tag plotting function.
    '''
    tags_data = {
        'agents': {
            time: timepoint['agents']
            for time, timepoint in data.items()
        },
        'config': environment_config,
    }
    return tags_data
