'''Generate figures for wcEcoli colony simulation

For usage information, run:
    python make_figures.py -h
'''
import argparse
from datetime import datetime
import json
import os
import sys
import subprocess
from typing import (
    Sequence, List, Dict, Tuple, Union, Iterable, Optional, Set, cast)

from matplotlib import rcParams  # type: ignore
import numpy as np
from vivarium.core.serialize import serialize_value
from vivarium.library.topology import get_in

from src.db import (
    add_connection_args,
    format_data_for_snapshots,
    format_data_for_tags,
    get_experiment_data,
)
from src.expression_survival import (
    plot_expression_survival,
    plot_expression_survival_dotplot,
)
from src.constants import OUT_DIR, FIELDS_PATH, BOUNDS_PATH
from src.types import RawData, EnvironmentConfig, SearchData, DataTuple
from src.total_mass import get_total_mass_plot
from src.environment_cross_sections import (
    get_enviro_sections_plot, SerializedField)
from src.phylogeny import plot_phylogeny
from src.process_expression_data import (
    raw_data_to_end_expression_table, VOLUME_KEY)
from src.ridgeline import get_ridgeline_plot
from src.plot_snapshots import plot_snapshots, plot_tags  # type: ignore
from src.centrality import get_survival_against_centrality_plot


# Colors from https://personal.sron.nl/~pault/#sec:qualitative
COLORS = {
    'blue': '#0077BB',
    'orange': '#EE7733',
    'cyan': '#33BBEE',
    'red': '#CC3311',
    'teal': '#009988',
    'magenta': '#EE3377'}
PUMP_PATH = (
    'periplasm', 'concentrations', 'TRANS-CPLX-201[m]')
BETA_LACTAMASE_PATH = (
    'periplasm', 'concentrations', 'EG10040-MONOMER[p]')
TAG_PATH_NAME_MAP = {
    ('bulk', 'EG10040-MONOMER[p]'): 'AmpC',
    ('bulk', 'TRANS-CPLX-201[m]'): 'AcrAB-TolC',
}
ENVIRONMENT_SECTION_FIELDS = ('GLC[p]',)
ENVIRONMENT_SECTION_TIMES: Tuple[int, ...] = (
    231, 6006, 11781, 17325, 23100)
AGENTS_TO_TRACE: Tuple[str, ...] = (
    '0_wcecoli01010101111',
    '0_wcecoli00100110',
    '0_wcecoli111111',
    '0_wcecoli01010100',
    '0_wcecoli0101011011',
    '0_wcecoli010110',
    '0_wcecoli1111010',
)
AGENTS_FOR_PHYLOGENY_TRACE = ('0_wcecoli101001101110',)
COLONY_MASS_PATH = ('mass',)
EXPRESSION_SURVIVAL_TIME_RANGE = (0.5, 1)
NUM_SNAPSHOTS = 5
FIG_OUT_DIR = os.path.join(OUT_DIR, 'figs')
FILE_EXTENSION = 'pdf'
ExperimentIdsType = Dict[
    str,
    Union[str, Tuple[str, ...], Dict[str, Tuple[str, ...]]],
]
EXPERIMENT_IDS: ExperimentIdsType = {
    'expression_distributions': (
        '20201119.150828', '20210112.185210', '20210125.152527'),
    'expression_heterogeneity': (
        '20201119.150828', '20210112.185210', '20210125.152527'),
    'enviro_heterogeneity': (
        '20201119.150828', '20210112.185210', '20210125.152527'),
    'enviro_section': (
        '20201119.150828', '20210112.185210', '20210125.152527'),
    'growth_basal': (
        '20201119.150828', '20210112.185210', '20210125.152527'),
    'growth_anaerobic': (
        '20201221.194828', '20210205.163802', '20210205.183800'),
    'growth': {
        'basal': (
            '20201119.150828', '20210112.185210', '20210125.152527'),
        'anaerobic': (
            '20201221.194828', '20210205.163802', '20210205.183800'),
    },
    'threshold_scan': {
        '0.01 mM': (
            '20201228.172246', '20210209.174715', '20210209.192621'),
        '0.02 mM': (
            '20201228.211700', '20210210.181817', '20210210.210112'),
        '0.03 mM': (
            '20201229.160649', '20210212.160943', '20210212.193824'),
        '0.04 mM': (
            '20201230.191552', '20210217.152358', '20210217.201244'),
        '0.05 mM': (
            '20210206.045834', '20210219.151929', '20210220.153658'),
    },
    'expression_survival': '20210329.155953',
    'expression_survival_dotplots': '20210329.155953',
    'death_snapshots': '20210329.155953',
    'death_snapshots_antibiotic': '20210329.155953',
    'centrality': '20210329.155953',
    'phylogeny': '20210329.155953',
}
FIGURE_NUMBER_NAME_MAP = {
    '3': {
        'A': 'enviro_heterogeneity',
        'B': 'enviro_section',
        'C': 'growth_basal',
        'D': 'growth_anaerobic',
        'E': 'growth',
        'F': 'expression_heterogeneity',
        'G': 'expression_distributions',
    },
    '5': {
        'A': 'threshold_scan',
        'B': 'death_snapshots',
        'C': 'centrality',
        'D': 'phylogeny',
        'E': 'expression_survival_dotplots',
        'F': 'expression_survival_dotplots',
        'G': 'expression_survival',
        'H': 'expression_survival',
        'I': 'expression_survival',
    },
    'X': {
        '1': 'death_snapshots_antibiotic',
    },
}
FIGURE_DESCRIPTIONS = {
    '3': {
        'A': 'snapshots of growing colony consuming glucose',
        'B': 'environment cross-sections showing glucose depletion',
        'C': 'snapshots in the basal condition',
        'D': 'snapshots in the anaerobic condition',
        'E': 'colony mass on basal and anaerobic media',
        'F': 'snapshots showing expression heterogeneity',
        'G': 'distributions of protein concentrations',
    },
    '5': {
        'A': 'parameter scan for tolerance threshold',
        'B': 'snapshot of final colony under nitrocefin',
        'C': 'box plot showing distances from center',
        'D': 'phylogenetic tree',
        'E': 'dotplot of final [AmpC] colored by survival',
        'F': 'dotplot of final [AcrAB-TolC] colored by survival',
        'G': 'final [AmpC] and plotted against [AcrAB-TolC]',
        'H': 'paths of dead agents through concentration space',
        'I': 'paths of a lineage of agents through concentration space',
    },
    'X': {
        '1': 'death_snapshots_antibiotic',
    },
}
METADATA_FILE = 'metadata.json'
STATS_FILE = 'stats.json'


def exec_shell(
        tokens: Sequence[str],
        timeout: int = 10,
        ) -> Tuple[str, str]:
    '''Execute a shell command and return the output.'''
    proc = subprocess.run(
        tokens,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
        env=None,
        universal_newlines=True,
        timeout=timeout)
    return proc.stdout.rstrip(), proc.stderr.rstrip()


def get_metadata() -> dict:
    '''Get information on which experiments and code were used.'''
    if os.environ.get('CI'):
        ci_url = '{}/{}/actions/runs/{}'.format(
            os.environ['GITHUB_SERVER_URL'],
            os.environ['GITHUB_REPOSITORY'],
            os.environ['GITHUB_RUN_ID'],
        )
        return {
            'git_hash': os.environ['GITHUB_SHA'],
            'git_branch': os.environ.get('GITHUB_BASE_REF'),
            'git_status': 'n/a because on CI',
            'time': datetime.utcnow().isoformat() + '+00:00',
            'python': sys.version.splitlines()[0],
            'experiment_ids': EXPERIMENT_IDS,
            'ci_url': ci_url,
        }
    return {
        'git_hash': exec_shell(['git', 'rev-parse', 'HEAD'])[0],
        'git_branch': exec_shell(
            ['git', 'symbolic-ref', '--short', 'HEAD'])[0],
        'git_status': exec_shell(
            ['git', 'status', '--porcelain'])[0].split('\n'),
        'time': datetime.utcnow().isoformat() + '+00:00',
        'python': sys.version.splitlines()[0],
        'experiment_ids': EXPERIMENT_IDS,
    }


def get_experiment_ids(
        id_obj: Union[str, Tuple[str, ...], List[str], Set[str], dict]
        ) -> List[str]:
    '''Get a flat list of all experiment IDs. May have duplicates.'''
    if isinstance(id_obj, str):
        return [id_obj]
    if isinstance(id_obj, (tuple, list, set)):
        ids_lst = []
        for elem in id_obj:
            ids_lst.extend(get_experiment_ids(elem))
        return ids_lst
    if isinstance(id_obj, dict):
        ids_lst = []
        for elem in id_obj.values():
            ids_lst.extend(get_experiment_ids(elem))
        return ids_lst
    return id_obj


def make_snapshots_figure(
        data: RawData,
        environment_config: EnvironmentConfig,
        name: str,
        fields: Sequence[str],
        agent_fill_color: Optional[str] = None,
        agent_alpha: float = 1,
        num_snapshots: int = NUM_SNAPSHOTS,
        snapshot_times: Optional[Tuple[float, ...]] = None,
        xlim: Tuple[float, float] = (0, 50),
        ylim: Tuple[float, float] = (0, 50)
        ) -> dict:
    '''Make a figure of snapshots.

    Args:
        data: The experiment data.
        environment_config: Environment parameters.
        name: Name of the output file (excluding file extension).
        fields: List of the names of fields to include.
        agent_fill_color: Fill color for agents.
        agent_alpha: Transparency for agents.
        num_snapshots: Number of snapshots.
        snapshot_times: Times to take snapshots at. If
            None, they are evenly spaced.
        xlim: Limits of x-axis.
        ylim: Limits of y-axis.

    Returns:
        Statistics.
    '''
    snapshots_data = format_data_for_snapshots(
        data, environment_config)
    if not fields:
        data = RawData({
            key: val
            for key, val in data.items() if key != 'fields'
        })
    plot_config = {
        'out_dir': FIG_OUT_DIR,
        'filename': '{}.{}'.format(name, FILE_EXTENSION),
        'include_fields': fields,
        'field_label_size': 54,
        'default_font_size': 54,
        'agent_fill_color': agent_fill_color,
        'dead_color': (0, 0, 0.79),  # gray in HSV
        'agent_alpha': agent_alpha,
        'n_snapshots': num_snapshots,
        'snapshot_times': snapshot_times,
        'scale_bar_length': 10,
        'scale_bar_color': 'white' if fields else 'black',
        'xlim': xlim,
        'ylim': ylim,
        'min_color': '#FFFFFF',
        'max_color': '#000000',
        'grid_color': 'white' if fields else '',
        'begin_gradient': 1,
    }
    stats = plot_snapshots(snapshots_data, plot_config)
    return stats


def make_expression_heterogeneity_fig(
        replicates_data: Iterable[DataTuple],
        _: SearchData,
        ) -> dict:
    '''Figure shows heterogeneous expression within wcEcoli agents.

    Create Figure 3F.
    '''
    for i, (data, enviro_config) in enumerate(replicates_data):
        tags_data = format_data_for_tags(data, enviro_config)
        tagged_molecules = list(TAG_PATH_NAME_MAP.keys())
        plot_config = {
            'out_dir': FIG_OUT_DIR,
            'tagged_molecules': tagged_molecules,
            'background_color': 'white',
            'filename': 'expression_heterogeneity_{}.{}'.format(
               i, FILE_EXTENSION),
            'tag_path_name_map': TAG_PATH_NAME_MAP,
            'tag_label_size': 54,
            'default_font_size': 48,
            'n_snapshots': NUM_SNAPSHOTS,
            'tag_colors': {
                tag: ('white', '#0000ff')
                for tag in tagged_molecules
            },
            'scale_bar_length': 10,
            'scale_bar_color': 'black',
            'xlim': (10, 40),
            'ylim': (10, 40),
        }
        plot_tags(tags_data, plot_config)
    return {}


def _calculate_distribution_stats(
        replicates_data: List[Tuple[Dict[str, Sequence[float]], str]],
        ) -> dict:
    stats = {}
    keys = replicates_data[0][0].keys()
    for replicate, _ in replicates_data:
        assert replicate.keys() == keys
    for key in keys:
        num_cells = 0
        key_replicates = []
        for replicate, _ in replicates_data:
            key_replicates.append(replicate[key])
            num_cells += len(replicate[key])
        key_replicates_array = np.array(key_replicates)  # type: ignore
        q1, q2, q3 = np.percentile(
            key_replicates_array,
            [25, 50, 75],
            axis=1,
        )
        stats[key] = (
            (key_replicates_array == 0).sum(),
            key_replicates_array.min(axis=1),  # type: ignore
            q1, q2, q3,
            key_replicates_array.max(axis=1),  # type: ignore
            num_cells,
        )
    return stats


def make_expression_distributions_fig(
        replicates_raw_data: Iterable[DataTuple],
        _search_data: SearchData
        ) -> dict:
    '''Figure shows the distributions of expression values.

    Create Figure 3G.
    '''
    replicates_data = []
    colors = ('#333333', '#777777', '#BBBBBB')
    for i, (raw_data, _) in enumerate(replicates_raw_data):
        color = colors[i]
        end_expression_table = raw_data_to_end_expression_table(
            raw_data,
            {val: key for key, val in TAG_PATH_NAME_MAP.items()})
        data = {
            key: end_expression_table[key].tolist()
            for key in end_expression_table.columns
            if key != VOLUME_KEY
        }
        replicates_data.append((data, color))
    stats = _calculate_distribution_stats(replicates_data)
    fig = get_ridgeline_plot(
        replicates_data,
        point_alpha=1,
        overlap=-0.1,
        horizontal_extra=0,
        num_bins=100,
        jitter=0,
        x_label='Protein Concentration (counts/fL)',
        y_label='Distribution Density',
        fontsize=14,
    )
    fig.savefig(
        os.path.join(
            FIG_OUT_DIR,
            'expression_distributions.{}'.format(FILE_EXTENSION),
        )
    )
    return stats


def make_growth_fig(
        raw_data: Dict[str, DataTuple],
        _: SearchData,
        ) -> dict:
    '''Make plot of colony mass of basal and anaerobic colonies.

    Create Figure 3E.
    '''
    data_dict = {
        'basal': [data for data, _ in raw_data['basal']],
        'anaerobic': [data for data, _ in raw_data['anaerobic']],
    }
    fig, stats = get_total_mass_plot(
        data_dict, list(COLORS.values()), fontsize=12)
    fig.savefig(os.path.join(
        FIG_OUT_DIR, 'growth.{}'.format(FILE_EXTENSION)))
    return stats


def make_threshold_scan_fig(
        data_and_configs: Dict[str, DataTuple],
        _: SearchData,
        ) -> dict:
    '''Plot colony mass curves with various antibiotic thresholds.

    Create Figure 5A.
    '''
    data_dict = dict({
        threshold: list(
            data for data, enviro_config in threshold_ids
        )
        for threshold, threshold_ids in data_and_configs.items()
    })
    some_data = list(data_dict.values())[0][0]
    vlines = (
        (
            max(some_data.keys()) / 2,
            0.51,
            'black',
            'Nitrocefin\nIntroduced',
        ),
    )
    fig, stats = get_total_mass_plot(
        data_dict, list(COLORS.values()), fontsize=12, vlines=vlines,
    )
    fig.savefig(os.path.join(
        FIG_OUT_DIR, 'threshold_scan.{}'.format(FILE_EXTENSION)))
    return stats


def make_expression_survival_fig(
        data_and_config: DataTuple,
        search_data: SearchData,
        ) -> dict:
    '''Make expression-survival figures.

    Create Figures 5G, 5H, and 5I.
    '''
    data, _ = data_and_config
    fig = plot_expression_survival(
        data, PUMP_PATH, BETA_LACTAMASE_PATH,
        'Final [AcrAB-TolC] (µM)',
        'Final [AmpC] (µM)',
        search_data['x_values'],
        search_data['y_values'],
        search_data['precision'],
        boundary_color=COLORS['magenta'],
        scaling=1e3,
        time_range=EXPRESSION_SURVIVAL_TIME_RANGE,
        fontsize=12,
    )
    fig.savefig(os.path.join(
        FIG_OUT_DIR, 'expression_survival.{}'.format(
            FILE_EXTENSION)
    ))
    plot_agents = set()
    for agent in AGENTS_FOR_PHYLOGENY_TRACE:
        for i in range(len('0_wcecoli') + 1, len(agent) + 1):
            plot_agents.add(agent[:i])
    fig = plot_expression_survival(
        data, PUMP_PATH, BETA_LACTAMASE_PATH,
        '[AcrAB-TolC] (µM)',
        '[AmpC] (µM)',
        search_data['x_values'],
        search_data['y_values'],
        search_data['precision'],
        boundary_color=COLORS['magenta'],
        scaling=1e3,
        time_range=EXPRESSION_SURVIVAL_TIME_RANGE,
        plot_agents=plot_agents,
        agents_for_phylogeny_trace=AGENTS_FOR_PHYLOGENY_TRACE,
        fontsize=12,
    )
    fig.savefig(os.path.join(
        FIG_OUT_DIR, 'expression_survival_lineage_traces.{}'.format(
            FILE_EXTENSION)
    ))
    fig = plot_expression_survival(
        data, PUMP_PATH, BETA_LACTAMASE_PATH,
        '[AcrAB-TolC] (µM)',
        '[AmpC] (µM)',
        search_data['x_values'],
        search_data['y_values'],
        search_data['precision'],
        boundary_color=COLORS['magenta'],
        scaling=1e3,
        time_range=EXPRESSION_SURVIVAL_TIME_RANGE,
        dead_trace_agents=AGENTS_TO_TRACE,
        plot_agents=AGENTS_TO_TRACE,
        fontsize=12,
    )
    fig.savefig(os.path.join(
        FIG_OUT_DIR, 'expression_survival_death_traces.{}'.format(
            FILE_EXTENSION)
    ))
    fig = plot_expression_survival(
        data, PUMP_PATH, BETA_LACTAMASE_PATH,
        'Final [AcrAB-TolC] (µM)',
        'Final [AmpC] (µM)',
        search_data['x_values'],
        search_data['y_values'],
        search_data['precision'],
        boundary_color=COLORS['magenta'],
        scaling=1e3,
        time_range=EXPRESSION_SURVIVAL_TIME_RANGE,
        label_agents=True,
        fontsize=12,
    )
    fig.savefig(os.path.join(
        FIG_OUT_DIR, 'expression_survival_labeled.{}'.format(
            FILE_EXTENSION)
    ))
    return {}


def make_expression_survival_dotplots(
        data_and_config: DataTuple,
        _search_data: SearchData,
        ) -> dict:
    '''Make dotplots of protein concentrations colored by survival.

    Create Figures 5E and 5F.
    '''
    data, _ = data_and_config
    stats = {}
    fig, stats['AcrAB-TolC'] = plot_expression_survival_dotplot(
        data, PUMP_PATH, 'Final [AcrAB-TolC] (µM)',
        scaling=1e3,
        time_range=EXPRESSION_SURVIVAL_TIME_RANGE,
        fontsize=12,
    )
    fig.savefig(os.path.join(
        FIG_OUT_DIR, 'expression_survival_pump.{}'.format(
            FILE_EXTENSION)))
    fig, stats['AmpC'] = plot_expression_survival_dotplot(
        data, BETA_LACTAMASE_PATH, 'Final [AmpC] (µM)',
        scaling=1e3,
        time_range=EXPRESSION_SURVIVAL_TIME_RANGE,
        fontsize=12,
    )
    fig.savefig(os.path.join(
        FIG_OUT_DIR,
        'expression_survival_beta_lactamase.{}'.format(
            FILE_EXTENSION)))
    return stats


def make_survival_centrality_fig(
        data_and_config: DataTuple,
        _search_data: SearchData,
        ) -> dict:
    '''Plot box plot figure of agent distances from center.

    Create Figure 5C.
    '''
    data, _ = data_and_config
    fig, stats = get_survival_against_centrality_plot(data)
    fig.savefig(os.path.join(
        FIG_OUT_DIR,
        'survival_centrality.{}'.format(FILE_EXTENSION)
    ))
    return stats


def make_environment_section(
        data_and_configs: Sequence[DataTuple],
        _search_data: SearchData,
        ) -> dict:
    '''Plot field concentrations in cross-section of final enviro.

    Create Figure 3B.
    '''
    t_final = max(data_and_configs[0][0].keys())
    fields_ts: List[Dict[float, Dict[str, SerializedField]]] = []
    section_times = [
        float(time) for time in ENVIRONMENT_SECTION_TIMES]
    for i, (replicate, _) in enumerate(data_and_configs):
        fields_ts.append(dict())
        for time in section_times:
            fields_ts[i][time] = {
                name: field
                for name, field in get_in(
                    replicate[time], FIELDS_PATH).items()
                if name in ENVIRONMENT_SECTION_FIELDS
            }
    bounds = get_in(data_and_configs[0][0][t_final], BOUNDS_PATH)
    fig, stats = get_enviro_sections_plot(
        fields_ts, bounds, section_location=0.5, fontsize=18)
    fig.savefig(
        os.path.join(FIG_OUT_DIR, 'enviro_section.{}'.format(
            FILE_EXTENSION)))
    return stats


def make_growth_basal_fig(
        replicates_data: Iterable[DataTuple],
        _: SearchData,
        ) -> dict:
    '''Create snapshots figure of colony on basal media.

    Create Figure 3C.
    '''
    stats = {}
    for i, (data, enviro_config) in enumerate(replicates_data):
        stats[i] = make_snapshots_figure(
            data, enviro_config, 'growth_basal_{}'.format(i), [])
    return stats


def make_growth_anaerobic_fig(
        replicates_data: Iterable[DataTuple],
        _: SearchData,
        ) -> dict:
    '''Create snapshots figure of colony on anaerobic media.

    Create Figure 3D.
    '''
    stats = {}
    for i, (data, enviro_config) in enumerate(replicates_data):
        stats[i] = make_snapshots_figure(
            data, enviro_config, 'growth_anaerobic_{}'.format(i), [])
    return stats


def make_phylogeny_plot(
        data_and_config: DataTuple,
        _search_data: SearchData,
        ) -> dict:
    '''Plot phylogenetic tree.

    Create Figure 5D.
    '''
    data, _ = data_and_config
    tree, df = plot_phylogeny(data, os.path.join(
        FIG_OUT_DIR, 'phylogeny.{}').format(FILE_EXTENSION),
        time_range=EXPRESSION_SURVIVAL_TIME_RANGE)
    tree.write(
        format=1,
        outfile=os.path.join(FIG_OUT_DIR, 'phylogeny.nw'),
    )
    df.to_csv(
        os.path.join(FIG_OUT_DIR, 'agent_survival.csv'),
        index=False,
    )
    return {}


def make_enviro_heterogeneity_fig(
        replicates_data: Iterable[DataTuple],
        _: SearchData,
        ) -> dict:
    '''Plot snapshots of colony consuming glucose.

    Create Figure 5A.
    '''
    stats = {}
    for i, (data, enviro_config) in enumerate(replicates_data):
        stats[i] = make_snapshots_figure(
            data, enviro_config, 'enviro_heterogeneity_{}'.format(i),
            ENVIRONMENT_SECTION_FIELDS, 'white',
        )
    return stats


def make_death_snapshots(
        data_and_config: DataTuple,
        _: SearchData,
        ) -> dict:
    '''Plot colony exposed to antibiotics with death coloration.

    Create Figure 5B.
    '''
    data, config = data_and_config
    return make_snapshots_figure(
        data, config, 'death_snapshots', [],
        agent_fill_color='green',
        snapshot_times=(max(data.keys()),),
        xlim=(5, 45),
        ylim=(5, 45),
    )

def make_death_snapshots_antibiotic(
        data_and_config: DataTuple,
        _: SearchData,
        ) -> dict:
    '''Plot colony consuming (or not) antibiotics.

    Create Figure X1. This figure is not used in the paper, but its
    statistics are used to find the change in nitrocefin concentration
    over time.
    '''
    data, config = data_and_config
    return make_snapshots_figure(
        data, config, 'death_snapshots_antibiotic',
        ['nitrocefin'],
        agent_fill_color='green',
        xlim=(5, 45),
        ylim=(5, 45),
    )


def create_data_dict(
        all_data: Dict[str, DataTuple],
        experiment_id_obj: Union[dict, str, Tuple[str, ...]],
        ) -> Union[dict, DataTuple, Tuple[DataTuple, ...]]:
    '''Create a dictionary of experiment simulation data.

    Note that we only support tuples of IDs, not tuples of dictionaries
    or tuples of tuples.

    Args:
        all_data (dict(str, tuple(RawData, dict))): Map from experiment
            ID to that experiment's simulation data and environment
            config.
        experiment_id_obj (Union(dict, str, tuple)): Object with one or
            more experiment IDs. The returned data object will have the
            same shape as this object. Dictionaries, strings, and tuples
            are supported.

    Returns: Data object of the same form as the experiment_id_obj.
    '''
    if isinstance(experiment_id_obj, str):
        return all_data[experiment_id_obj]
    if isinstance(experiment_id_obj, tuple):
        to_return = []
        for elem in experiment_id_obj:
            assert isinstance(elem, str)
            data_tuple = create_data_dict(all_data, elem)
            to_return.append(cast(DataTuple, data_tuple))
        return tuple(to_return)
    if isinstance(experiment_id_obj, dict):
        return dict({
            key: create_data_dict(all_data, value)
            for key, value in experiment_id_obj.items()
        })
    raise ValueError(
        'experiment_id_obj %s of unsupported type %s' %
        (experiment_id_obj, type(experiment_id_obj)),
    )


FIGURE_FUNCTION_MAP = {
    'expression_distributions': make_expression_distributions_fig,
    'expression_heterogeneity': make_expression_heterogeneity_fig,
    'growth_basal': make_growth_basal_fig,
    'growth_anaerobic': make_growth_anaerobic_fig,
    'growth': make_growth_fig,
    'enviro_heterogeneity': make_enviro_heterogeneity_fig,
    'enviro_section': make_environment_section,
    'threshold_scan': make_threshold_scan_fig,
    'expression_survival': make_expression_survival_fig,
    'expression_survival_dotplots': make_expression_survival_dotplots,
    'phylogeny': make_phylogeny_plot,
    'death_snapshots': make_death_snapshots,
    'death_snapshots_antibiotic': make_death_snapshots_antibiotic,
    'centrality': make_survival_centrality_fig,
}


def main() -> None:
    '''Generate all figures.'''
    rcParams['font.sans-serif'] = ['Arial']
    rcParams['font.family'] = ['sans-serif']
    if not os.path.exists(FIG_OUT_DIR):
        os.makedirs(FIG_OUT_DIR)
    with open(os.path.join(FIG_OUT_DIR, METADATA_FILE), 'w') as f:
        json.dump(get_metadata(), f, indent=4)
    stats: dict = {}
    parser = argparse.ArgumentParser(
        description=(
            'Generate selected figures and associated stats from '
            'simulation data.'
        ),
    )
    add_connection_args(parser)
    parser.add_argument(
        'search_data', type=str, help='Path to boundary search data.')
    parser.add_argument(
        '--data_path',
        default='',
        help='Folder of JSON files to read data from instead of Mongo',
    )
    for figure, d in FIGURE_NUMBER_NAME_MAP.items():
        for panel in d:
            parser.add_argument(
                '--{}{}'.format(figure, panel),
                action='store_true',
                help='Generate figure & stats for fig {}{}: {}'.format(
                    figure, panel, FIGURE_DESCRIPTIONS[figure][panel]),
            )
    parser.add_argument(
        '--all',
        action='store_true',
        help='Generate all figures and stats.',
    )
    args = parser.parse_args()
    args_dict = vars(args)

    with open(args.search_data, 'r') as f:
        search_data = json.load(f)

    data_cache = {}
    generated_figures = set()
    stats = {}

    for fig, fig_dict in FIGURE_NUMBER_NAME_MAP.items():
        for panel, fig_name in fig_dict.items():
            if not (args_dict['{}{}'.format(fig, panel)] or args.all):
                continue
            if fig_name in generated_figures:
                continue
            generated_figures.add(fig_name)
            experiment_ids = EXPERIMENT_IDS[fig_name]
            for experiment_id in get_experiment_ids(experiment_ids):
                if experiment_id in data_cache:
                    continue
                data_cache[experiment_id] = get_experiment_data(
                    args, experiment_id)
            data = create_data_dict(data_cache, experiment_ids)
            func = FIGURE_FUNCTION_MAP[fig_name]
            stats[fig_name] = func(data, search_data)  # type: ignore

    with open(os.path.join(FIG_OUT_DIR, STATS_FILE), 'w') as f:
        json.dump(serialize_value(stats), f, indent=4)


if __name__ == '__main__':
    main()
