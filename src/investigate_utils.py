'''Utilities for investigation scripts.
'''

import argparse
import copy
import os

from typing import Dict, Tuple
from vivarium.core.experiment import get_in, assoc_path
from vivarium_cell.analysis.analyze import Analyzer

from src.types import RawData
from src.constants import OUT_DIR


PATH_TO_AGENTS = ('agents',)
PATH_TO_DEAD = ('boundary', 'dead')


def filter_raw_data_by_time(
        raw_data: RawData,
        time_range: Tuple[float, float]) -> RawData:
    '''Filter raw simulation data to the timepoints within a range

    Args:
        raw_data: Raw simulation data.
        time_range: Tuple of range endpoints. Each endpoint is a float
            between 0 and 1 (inclusive) that denotes a fraction of the
            total simulation time.
    Returns:
        A subset of the key-value pairs in ``raw_data``. Includes only
        those timepoints between the ``time_range`` endpoints
        (inclusive).
    '''
    floor, ceil = time_range
    end = max(raw_data.keys())
    filtered = RawData({
        time: time_data
        for time, time_data in raw_data.items()
        if floor * end <= time <= ceil * end
    })
    return filtered


def split_raw_data_by_survival(
        raw_data: RawData) -> Tuple[RawData, RawData]:
    '''Segment raw data into data for agents that die and that survive

    Args:
        raw_data: Raw simulation data
    Returns:
        Tuple of 2 raw data dictionaries. The first contains all agents
        that survive until division. The second contains all agents that
        die before dividing.
    '''
    # Establish which agents die
    agents_die = set()
    for time_data in raw_data.values():
        agents_data = get_in(time_data, PATH_TO_AGENTS)
        for agent, agent_data in agents_data.items():
            dead = get_in(agent_data, PATH_TO_DEAD, False)
            if dead:
                agents_die.add(agent)

    # Split the data
    survive_data = RawData(dict())
    for time in raw_data:
        agents_path = (time,) + PATH_TO_AGENTS
        assoc_path(survive_data, agents_path, dict())
    die_data = copy.deepcopy(survive_data)

    for time, time_data in raw_data.items():
        agents_data = get_in(time_data, PATH_TO_AGENTS)
        for agent, agent_data in agents_data.items():
            dest = die_data if agent in agents_die else survive_data
            agent_path = (time,) + PATH_TO_AGENTS + (agent,)
            assoc_path(dest, agent_path, agent_data)

    return survive_data, die_data


def add_experiment_id_arg(
        parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    '''Add experiment ID argument to parser.

    Args:
        parser: Argument parser.
    Returns:
        The argument parser with the experiment ID argument added.
    '''
    parser.add_argument(
        'experiment_id',
        type=str,
        help='ID of experiment to analyze.'
    )
    return parser


def create_dir(path: str) -> None:
    '''Create directory if it doesn't already exist.

    Also creates any intermediate directories.

    Args:
        path: Path to directory to create.
    '''
    if not os.path.exists(path):
        os.makedirs(path)


def parse_args_retrieve_data_by_id() -> Tuple[RawData, Dict, str]:
    '''Handle argument parsing and retrieve data by experiment ID

    For simple scripts that don't need to add any more connection
    arguments, this function handles adding the experiment ID, parsing
    the arguments, creating the output directory, and retrieving the
    data.

    Returns:
        Tuple of the raw experiment data, environment configuration, and
        output directory path.
    '''
    parser = argparse.ArgumentParser()
    Analyzer.add_connection_args(parser)
    add_experiment_id_arg(parser)
    args = parser.parse_args()

    out_dir = os.path.join(OUT_DIR, args.experiment_id)
    create_dir(out_dir)

    data, environment_config = Analyzer.get_data(
        args, args.experiment_id)
    return data, environment_config, out_dir
