import argparse
import copy
import json
import os
from typing import Union, Optional, cast

from vivarium.core.emitter import (
    get_local_client,
    data_from_database,
)
from vivarium.core.serialize import deserialize_value
from vivarium.library.units import remove_units

from src.types import RawData, EnvironmentConfig, DataTuple


class SplitExperimentSpec:

    def __init__(self, spec: dict):
        self._spec = copy.deepcopy(spec)

    def items(self):
        return self._spec.items()

    def key(self):
        return tuple(self.items())

    def to_dict(self):
        return copy.deepcopy(self._spec)

    def __hash__(self):
        return hash(self.key())

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.key() == other.key()


def get_experiment_data(
        args: argparse.Namespace,
        experiment_id_spec: Union[str, SplitExperimentSpec],
        ) -> DataTuple:
    '''Get simulation data for an experiment.

    If ``args.data_path`` is set, retrieve the experiment data from a
    JSON file named ``<experiment_id>.json`` under ``args.data_path``.
    Otherwise, retrieve the data from MongoDB.

    Args:
        args: Parsed CLI args.
        experiment_id_spec: Specification for the ID(s) of the
            experiment(s) to retrieve data for. If a string, it should
            be the experiment ID. If a dictionary, it should be a
            mapping from start times to experiment IDs (strings).

    Returns: Tuple of simulation data and environment config. If
        ``experiment_id_spec`` was a dictionary, only the first
        environment config will be returned, and the simulation data
        will contain data from all the experiments, with each time
        increased by the experiment's start time.

    Raises:
        AssertionError: If there are overlapping times between
            experiments or if there is no environment config.
    '''
    if isinstance(experiment_id_spec, str):
        experiment_ids = SplitExperimentSpec({
            0.0: experiment_id_spec,
        })
    else:
        experiment_ids = cast(SplitExperimentSpec, experiment_id_spec)

    environment_config: Optional[EnvironmentConfig] = None
    combined_data: dict = {}
    for start_time, experiment_id in experiment_ids.items():
        path = os.path.join(
            args.data_path, '{}.json'.format(experiment_id))
        if args.data_path and os.path.exists(path):
            with open(path, 'r') as f:
                loaded_file = json.load(f)
                data = RawData({
                    float(time) + start_time: value
                    for time, value in loaded_file['data'].items()
                })
                if environment_config is None:
                    environment_config = EnvironmentConfig(
                        loaded_file['environment_config'])
                assert not combined_data.keys() & data.keys()
                combined_data.update(data)
        else:
            client = get_local_client(
                args.host, args.port, args.database_name)
            data, _ = data_from_database(
                experiment_id, client,
                #filters={'data.time': {'$mod': [10, 0]}}
            )
            data = remove_units(deserialize_value(data))
            if environment_config is None:
                environment_config = data[min(data)]['dimensions']
            data = RawData({
                time + start_time: timepoint
                for time, timepoint in data.items()
            })
            assert not combined_data.keys() & data.keys()
            combined_data.update(data)
    assert environment_config is not None
    return RawData(combined_data), environment_config


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
