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
    SplitExperimentSpec,
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
from src.types import (
    RawData,
    EnvironmentConfig,
    SearchData,
    DataTuple,
)
from src.total_mass import get_total_mass_plot
from src.environment_cross_sections import (
    get_enviro_sections_plot, SerializedField)
from src.phylogeny import plot_phylogeny
from src.process_expression_data import (
    raw_data_to_end_expression_table, VOLUME_KEY)
from src.ridgeline import get_ridgeline_plot
from src.plot_snapshots import plot_snapshots, plot_tags  # type: ignore
from src.centrality import get_survival_against_centrality_plot
from src.timeseries import get_timeseries_plot
from src.vivarium_wcecoli_comparison import (
    get_proteome_comparison_plot,
    get_submass_comparison_plot,
    get_mass_fraction_comparison_plot,
)


# Colors from https://personal.sron.nl/~pault/#sec:qualitative
COLORS = {
    'blue': '#0077BB',
    'orange': '#EE7733',
    'cyan': '#33BBEE',
    'red': '#CC3311',
    'teal': '#009988',
    'magenta': '#EE3377'
}
COLOR_LIST = list(COLORS.values())
PUMP_PATH = (
    'periplasm', 'concentrations', 'TRANS-CPLX-201[m]')
BETA_LACTAMASE_PATH = (
    'periplasm', 'concentrations', 'EG10040-MONOMER[p]')
TAG_PATH_NAME_MAP = {
    ('bulk', 'EG10040-MONOMER[p]'): 'AmpC',
    ('bulk', 'TRANS-CPLX-201[m]'): 'AcrAB-TolC',
}
EQUILIBRATION_TIME_RANGE = (25, 35)
EQUILIBRATION_AGENT = '0'
EQUILIBRATION_HORIZONTAL_LINES = {
    'tetracycline': (
        ('external', 0.0025, COLOR_LIST[1]),
        # Source: Thanassi et al. (1995)
        ('expected', 10e-3, COLOR_LIST[2]),
    ),
    'ampicillin': (
        ('external', 0.0057, COLOR_LIST[1]),
        # Source: Kojima and Nikaido (2013)
        ('expected', 1.7e-3, COLOR_LIST[2]),
    ),
}
GROWTH_RATE_AGENT = '0'
RNA_PROTEIN_NAME_PATHS_MAP = {
    'AmpC': (
        (('bulk', 'EG10040-MONOMER[p]'), 1),
    ),
    'ampC': (
        (('bulk', 'EG10040_RNA[c]'), 1),
    ),
    'AcrA': (
        (('bulk', 'EG11703-MONOMER[p]'), 1),
        (('bulk', 'TRANS-CPLX-201[m]'), 6),
        (('bulk', 'CPLX0-3932[i]'), 6),
    ),
    'acrA': (
        (('bulk', 'EG11703_RNA[c]'), 1),
    ),
    'MarR': (
        (('bulk', 'PD00364[c]'), 1),
        (('bulk', 'CPLX0-7710[c]'), 2),
    ),
    'marR': (
        (('bulk', 'EG11435_RNA[c]'), 1),
    ),
    'MarA': (
        (('bulk', 'PD00365[c]'), 1),
    ),
    'marA': (
        (('bulk', 'EG11434_RNA[c]'), 1),
    ),
}
RNA_PROTEIN_LAYOUT = [
    ['acrA', 'ampC', 'marA', 'marR'],
    ['AcrA', 'AmpC', 'MarA', 'MarR'],
]
RNA_PROTEIN_AGENT = '00'
TETRACYCLINE_ACTIVITY_TIMESERIES_NAME_PATHS_MAP = {
    'Inhibited Active Ribosomes': (
        (('listeners', 'unique_counts', 'active_ribosome_tetracycline'), 1),
    ),
    'Inhibited 30S Subunits': (
        (('bulk', 'CPLX0-3953-tetracycline[c]'), 1),
    ),
    'Uninhibited 30S Subunits': (
        (('bulk', 'CPLX0-3953[c]'), 1),
    ),
    'Uninhibited Active Ribosomes': (
        (('listeners', 'unique_counts', 'active_ribosome'), 1),
    ),
}
TETRACYCLINE_TRANSPORT_TIMESERIES_NAME_PATHS_MAP = {
    'Cytoplasmic Tetracycline': (
        (('cytoplasm', 'concentrations', 'tetracycline'), 1),
    ),
    'External Tetracycline': (
        (('boundary', 'external', 'tetracycline'), 1),
    ),
    'OmpF': (
        (('bulk', 'CPLX0-7534[o]'), 1),
    ),
    'Outer Membrane': (
        (('kinetic_parameters', 'outer_tetracycline_permeability'), 1),
    ),
    'AcrAB-TolC': (
        (('bulk', 'TRANS-CPLX-201[m]'), 1),
    ),
}
TETRACYCLINE_TRANSPORT_TIMESERIES_Y_LABELS = {
    'Cytoplasmic Tetracycline': 'concentration (mM)',
    'External Tetracycline': 'concentration (mM)',
    'OmpF': 'counts',
    'Outer Membrane': 'permeability (cm/s)',
    'AcrAB-TolC': 'concentration (mM)',
}
TETRACYCLINE_TRANSPORT_TIMESERIES_LAYOUT = [
    ['External Tetracycline', 'Cytoplasmic Tetracycline', None],
    ['OmpF', 'Outer Membrane', 'AcrAB-TolC'],
]
TETRACYCLINE_ACTIVITY_TIMESERIES_LAYOUT = [
    ['Inhibited Active Ribosomes', 'Uninhibited Active Ribosomes'],
    ['Inhibited 30S Subunits', 'Uninhibited 30S Subunits'],
]
TETRACYCLINE_TRANSPORT_TIMESERIES_AGENT = '0'
TETRACYCLINE_ACTIVITY_TIMESERIES_AGENT = '0'
ENVIRONMENT_SECTION_FIELDS = ('GLC[p]',)
ENVIRONMENT_SECTION_TIMES: Tuple[int, ...] = (
    0, 2890, 5780, 8660, 11550)
AGENTS_TO_TRACE: Tuple[str, ...] = (
    '010',
)
AGENTS_FOR_PHYLOGENY_TRACE = ('010',)
COLONY_MASS_PATH = ('mass',)
EXPRESSION_SURVIVAL_TIME_RANGE = (0.5, 1)
NUM_SNAPSHOTS = 5
FIG_OUT_DIR = os.path.join(OUT_DIR, 'figs')
FILE_EXTENSION = 'pdf'
ExperimentIdsType = Dict[
    str,
    Union[str, tuple, dict],
]
EXPERIMENT_IDS: ExperimentIdsType = {
    'rna_protein_timeseries': (
        '709652e8-03b7-11ed-acfe-2f08daca2550'),
    'equilibration_tetracycline': (
        '0c44f61c-0477-11ed-acfe-2f08daca2550'),
    'growth_rate': (
        '447e44ba-055e-11ed-a7c0-87492988b953'),
    'tetracycline_activity_timeseries': (
        '98e8a2c0-0621-11ed-80f5-fbbd58ffc717'),
    'tetracycline_transport_timeseries': (
        '0c44f61c-0477-11ed-acfe-2f08daca2550'),
    'equilibration_ampicillin': (
        'bb92b112-049b-11ed-acfe-2f08daca2550'),
    'expression_distributions': (
        '1b2a3ca2-fd4f-11ec-ae52-8da6b112d368',),
    'expression_heterogeneity': (
        SplitExperimentSpec({
            0: '1b2a3ca2-fd4f-11ec-ae52-8da6b112d368',
            11551: '2022-08-11_22-31-14_068632+0000',
        }),
    ),
    'enviro_heterogeneity': (
        SplitExperimentSpec({
            0: '1b2a3ca2-fd4f-11ec-ae52-8da6b112d368',
            11551: '2022-08-11_22-31-14_068632+0000',
        }),
    ),
    'enviro_section': (
        SplitExperimentSpec({
            0: '1b2a3ca2-fd4f-11ec-ae52-8da6b112d368',
            11551: '2022-08-11_22-31-14_068632+0000',
        }),
    ),
    'growth_basal': (
        SplitExperimentSpec({
            0: '1b2a3ca2-fd4f-11ec-ae52-8da6b112d368',
            11551: '2022-08-11_22-31-14_068632+0000',
        }),
    ),
    'growth_anaerobic': (
        '4de8a748-bddd-11ec-af73-a566a8a29bc0',),
    'growth': {
        'basal': (
            SplitExperimentSpec({
                0: '1b2a3ca2-fd4f-11ec-ae52-8da6b112d368',
                11551: '2022-08-11_22-31-14_068632+0000',
            }),
        ),
        'anaerobic': (
            'fc32e70c-eb49-11ec-8f84-bb7045ed26ca',),
    },
    'threshold_scan': {
        '0.01 mM': (
            '27c2aa52-c12d-11ec-a3cd-d588635fd3e2',),
    },
    'expression_survival': '27c2aa52-c12d-11ec-a3cd-d588635fd3e2',
    'expression_survival_dotplots': '27c2aa52-c12d-11ec-a3cd-d588635fd3e2',
    'death_snapshots': '27c2aa52-c12d-11ec-a3cd-d588635fd3e2',
    'death_snapshots_antibiotic': '27c2aa52-c12d-11ec-a3cd-d588635fd3e2',
    'centrality': '27c2aa52-c12d-11ec-a3cd-d588635fd3e2',
    'phylogeny': '27c2aa52-c12d-11ec-a3cd-d588635fd3e2',
    'proteome_comparison': {
        'vivarium-ecoli': (
            '709652e8-03b7-11ed-acfe-2f08daca2550',
            'ffe868bc-201a-11ed-881c-a91fcb0e2640',
            'd5956e52-202e-11ed-881c-a91fcb0e2640',
            '30aa30b0-224c-11ed-881c-a91fcb0e2640',
            'c096dde6-225f-11ed-881c-a91fcb0e2640',
        ),
        'wcecoli': (
            'wcecoli_sim_seed_0',
            'wcecoli_sim_seed_1',
            'wcecoli_sim_seed_2',
            'wcecoli_sim_seed_3',
            'wcecoli_sim_seed_4',
        ),
    },
    'submass_comparison': {
        'vivarium-ecoli': (
            '709652e8-03b7-11ed-acfe-2f08daca2550',
            'ffe868bc-201a-11ed-881c-a91fcb0e2640',
            'd5956e52-202e-11ed-881c-a91fcb0e2640',
            '30aa30b0-224c-11ed-881c-a91fcb0e2640',
            'c096dde6-225f-11ed-881c-a91fcb0e2640',
        ),
        'wcecoli': (
            'wcecoli_sim_seed_0',
            'wcecoli_sim_seed_1',
            'wcecoli_sim_seed_2',
            'wcecoli_sim_seed_3',
            'wcecoli_sim_seed_4',
        ),
    },
    'mass_fraction_comparison': {
        'vivarium-ecoli': (
            '709652e8-03b7-11ed-acfe-2f08daca2550',
            'ffe868bc-201a-11ed-881c-a91fcb0e2640',
            'd5956e52-202e-11ed-881c-a91fcb0e2640',
            '30aa30b0-224c-11ed-881c-a91fcb0e2640',
            'c096dde6-225f-11ed-881c-a91fcb0e2640',
        ),
        'wcecoli': (
            'wcecoli_sim_seed_0',
            'wcecoli_sim_seed_1',
            'wcecoli_sim_seed_2',
            'wcecoli_sim_seed_3',
            'wcecoli_sim_seed_4',
        ),
    },
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
        '2': 'rna_protein_timeseries',
        '3': 'equilibration_tetracycline',
        '4': 'equilibration_ampicillin',
        '5': 'growth_rate',
        '6': 'tetracycline_activity_timeseries',
        '7': 'tetracycline_transport_timeseries',
        '8': 'proteome_comparison',
        '9': 'submass_comparison',
        '10': 'mass_fraction_comparison',
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
        '2': 'timeseries of selected RNA and protein counts',
        '3': 'timeseries of tetracycline reaching equilibrium',
        '4': 'timeseries of ampicillin reaching equilibrium',
        '5': 'timeseries of instantaneous growth rate',
        '6': 'timeseries of variables related to tetracycline activity',
        '7': 'timeseries of variables related to tetracycline transport',
        '8': 'comparison of wcEcoli and vivarium-ecoli proteomes',
        '9': 'comparison of wcEcoli and vivarium-ecoli submass growth',
        '10': 'comparison of wcEcoli and vivarium-ecoli mass fractions',
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


def _prepare_for_json(experiment_ids):
    if isinstance(experiment_ids, SplitExperimentSpec):
        return _prepare_for_json(experiment_ids.to_dict())
    if isinstance(experiment_ids, tuple):
        return tuple(_prepare_for_json(elem) for elem in experiment_ids)
    if isinstance(experiment_ids, str):
        return experiment_ids
    if isinstance(experiment_ids, dict):
        return {
            str(key): _prepare_for_json(value)
            for key, value in experiment_ids.items()
        }
    raise ValueError(
        f'Cannot prepare value of type {type(experiment_ids)} '
        f'for JSON: {experiment_ids}.')


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
            'experiment_ids': _prepare_for_json(EXPERIMENT_IDS),
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
        'experiment_ids': _prepare_for_json(EXPERIMENT_IDS),
    }


def get_experiment_ids(
        id_obj: Union[str, Tuple[str, ...], List[str], Set[str], dict]
        ) -> List[Union[str, SplitExperimentSpec]]:
    '''Get a flat list of all experiment IDs. May have duplicates.'''
    if isinstance(id_obj, str):
        return [id_obj]
    if isinstance(id_obj, (tuple, list, set)):
        ids_lst = []
        for elem in id_obj:
            ids_lst.extend(get_experiment_ids(elem))
        return ids_lst
    if isinstance(id_obj, SplitExperimentSpec):
        return [id_obj]
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
        for i in range(1, len(agent) + 1):
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


def make_rna_protein_timeseries(
        data_and_config: DataTuple,
        _: SearchData,
        ) -> dict:
    '''Plot RNA and protein expression timeseries.

    Create Figure X2. This figure is not used in the paper.
    '''
    data, _config = data_and_config
    fig = get_timeseries_plot(
        data,
        RNA_PROTEIN_NAME_PATHS_MAP,
        {name: 'count' for name in RNA_PROTEIN_NAME_PATHS_MAP},
        RNA_PROTEIN_LAYOUT,
        agent_path=('agents', RNA_PROTEIN_AGENT),
    )
    out_path = os.path.join(
        FIG_OUT_DIR, f'rna_protein_timeseries.{FILE_EXTENSION}')
    fig.savefig(out_path)
    return {}


def make_tetracycline_equilibration_timeseries(
        data_and_config: DataTuple,
        _: SearchData,
        ) -> dict:
    '''Plot timeseries of tetracycline equilibrating with environment.

    Create Figure X3. This figure is not used in the paper.
    '''
    data, _config = data_and_config
    min_time, max_time = EQUILIBRATION_TIME_RANGE
    data = RawData({
        time: timepoint
        for time, timepoint in data.items()
        if min_time <= time <= max_time
    })
    fig = get_timeseries_plot(
        data,
        {
            'tetracycline': (
                (('cytoplasm', 'concentrations', 'tetracycline'), 1),
            ),
        },
        {'tetracycline': 'concentration (mM)'},
        [['tetracycline']],
        col_width=8,
        row_height=4,
        agent_path=('agents', EQUILIBRATION_AGENT),
        horizontal_lines=EQUILIBRATION_HORIZONTAL_LINES,
    )
    out_path = os.path.join(
        FIG_OUT_DIR, f'equilibration_tetracycline.{FILE_EXTENSION}')
    fig.savefig(out_path)
    return {}


def make_ampicillin_equilibration_timeseries(
        data_and_config: DataTuple,
        _: SearchData,
        ) -> dict:
    '''Plot timeseries of ampicillin equilibrating with environment.

    Create Figure X4. This figure is not used in the paper.
    '''
    data, _config = data_and_config
    min_time, max_time = EQUILIBRATION_TIME_RANGE
    data = RawData({
        time: timepoint
        for time, timepoint in data.items()
        if min_time <= time <= max_time
    })
    fig = get_timeseries_plot(
        data,
        {
            'ampicillin': (
                (('periplasm', 'concentrations', 'ampicillin'), 1),
            ),
        },
        {'ampicillin': 'concentration (mM)'},
        [['ampicillin']],
        col_width=8,
        row_height=4,
        agent_path=('agents', EQUILIBRATION_AGENT),
        horizontal_lines=EQUILIBRATION_HORIZONTAL_LINES,
    )
    out_path = os.path.join(
        FIG_OUT_DIR, f'equilibration_ampicillin.{FILE_EXTENSION}')
    fig.savefig(out_path)
    return {}


def make_growth_rate_timeseries(
        data_and_config: DataTuple,
        _: SearchData,
        ) -> dict:
    '''Plot timeseries of growth rate under tetracycline.

    Create Figure X5. This figure is not used in the paper.
    '''
    data, _config = data_and_config
    fig = get_timeseries_plot(
        data,
        {
            'growth rate': (
                (('listeners', 'mass', 'instantaniousGrowthRate'), 1),
            ),
        },
        {'growth rate': 'fg/s'},
        [['growth rate']],
        col_width=8,
        row_height=4,
        agent_path=('agents', GROWTH_RATE_AGENT),
        min_value=0,
    )
    out_path = os.path.join(
        FIG_OUT_DIR, f'growth_rate.{FILE_EXTENSION}')
    fig.savefig(out_path)
    return {}


def make_tetracycline_activity_timeseries(
        data_and_config: DataTuple,
        _: SearchData,
        ) -> dict:
    data, _config = data_and_config
    fig = get_timeseries_plot(
        data,
        TETRACYCLINE_ACTIVITY_TIMESERIES_NAME_PATHS_MAP,
        {
            var: 'counts' for var in
            TETRACYCLINE_ACTIVITY_TIMESERIES_NAME_PATHS_MAP
        },
        TETRACYCLINE_ACTIVITY_TIMESERIES_LAYOUT,
        agent_path=('agents', TETRACYCLINE_ACTIVITY_TIMESERIES_AGENT),
    )
    out_path = os.path.join(
        FIG_OUT_DIR, f'tetracycline_activity_timeseries.{FILE_EXTENSION}')
    fig.savefig(out_path)
    return {}


def make_tetracycline_transport_timeseries(
        data_and_config: DataTuple,
        _: SearchData,
        ) -> dict:
    data, _config = data_and_config
    fig = get_timeseries_plot(
        data,
        TETRACYCLINE_TRANSPORT_TIMESERIES_NAME_PATHS_MAP,
        TETRACYCLINE_TRANSPORT_TIMESERIES_Y_LABELS,
        TETRACYCLINE_TRANSPORT_TIMESERIES_LAYOUT,
        agent_path=('agents', TETRACYCLINE_TRANSPORT_TIMESERIES_AGENT),
    )
    out_path = os.path.join(
        FIG_OUT_DIR, f'tetracycline_transport_timeseries.{FILE_EXTENSION}')
    fig.savefig(out_path)
    return {}


def make_proteome_comparison(
        data_and_configs: Dict[str, Iterable[DataTuple]],
        _: SearchData,
        ) -> dict:
    data = {
        key: [datatuple[0] for datatuple in datatuples]
        for key, datatuples in data_and_configs.items()
    }
    fig = get_proteome_comparison_plot(
        data,
        vivarium_agent='0',
    )
    out_path = os.path.join(
        FIG_OUT_DIR, f'proteome_comparison.{FILE_EXTENSION}')
    fig.savefig(out_path)
    return {}


def make_submass_comparison(
        data_and_configs: Dict[str, Iterable[DataTuple]],
        _: SearchData,
        ) -> dict:
    data = {
        key: [datatuple[0] for datatuple in datatuples]
        for key, datatuples in data_and_configs.items()
    }
    fig = get_submass_comparison_plot(
        data,
        COLOR_LIST,
        vivarium_agent='0',
    )
    out_path = os.path.join(
        FIG_OUT_DIR, f'submass_comparison.{FILE_EXTENSION}')
    fig.savefig(out_path)
    return {}


def make_mass_fraction_comparison(
        data_and_configs: Dict[str, Iterable[DataTuple]],
        _: SearchData,
        ) -> dict:
    data = {
        key: [datatuple[0] for datatuple in datatuples]
        for key, datatuples in data_and_configs.items()
    }
    fig = get_mass_fraction_comparison_plot(
        data,
        COLOR_LIST,
        vivarium_agent='0',
    )
    out_path = os.path.join(
        FIG_OUT_DIR, f'mass_fraction_comparison.{FILE_EXTENSION}')
    fig.savefig(out_path)
    return {}


def create_data_dict(
        all_data: Dict[Union[str, SplitExperimentSpec], DataTuple],
        experiment_id_obj: Union[
            dict, str, Tuple[str, ...],
            Tuple[SplitExperimentSpec, ...], SplitExperimentSpec],
        ) -> Union[dict, DataTuple, Tuple[DataTuple, ...]]:
    '''Create a dictionary of experiment simulation data.

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
    if isinstance(experiment_id_obj, (str, SplitExperimentSpec)):
        return all_data[experiment_id_obj]
    if isinstance(experiment_id_obj, tuple):
        to_return = []
        for elem in experiment_id_obj:
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
    'rna_protein_timeseries': make_rna_protein_timeseries,
    'equilibration_tetracycline': (
        make_tetracycline_equilibration_timeseries),
    'equilibration_ampicillin': (
        make_ampicillin_equilibration_timeseries),
    'growth_rate': make_growth_rate_timeseries,
    'tetracycline_activity_timeseries': (
        make_tetracycline_activity_timeseries),
    'tetracycline_transport_timeseries': (
        make_tetracycline_transport_timeseries),
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
    'proteome_comparison': make_proteome_comparison,
    'submass_comparison': make_submass_comparison,
    'mass_fraction_comparison': make_mass_fraction_comparison,
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
