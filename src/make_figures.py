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

from matplotlib import colors as mcolors
from vivarium.core.experiment import get_in
from vivarium_cell.analysis.analyze import Analyzer

from src.expression_survival import (
    plot_expression_survival,
    plot_expression_survival_dotplot,
)
from src.constants import OUT_DIR, FIELDS_PATH, BOUNDS_PATH
from src.total_mass import get_total_mass_plot
from src.environment_cross_sections import get_enviro_sections_plot
from src.phylogeny import plot_phylogeny
from src.process_expression_data import (
    raw_data_to_end_expression_table, VOLUME_KEY)
from src.ridgeline import get_ridgeline_plot
from src.expression_survival_scan import (
    load_scan_data, plot_expression_survival_scan)
from src.plot_snapshots import plot_snapshots, plot_tags


# Colors from https://personal.sron.nl/~pault/#sec:qualitative
COLORS = {
    'blue': '#0077BB',
    'orange': '#EE7733',
    'cyan': '#33BBEE',
    'red': '#CC3311',
    'teal': '#009988',
    'magenta': '#EE3377'}
PUMP_PATH = (
    'boundary', 'bulk_molecule_concentrations', 'TRANS-CPLX-201[s]')
BETA_LACTAMASE_PATH = (
    'boundary', 'bulk_molecule_concentrations', 'EG10040-MONOMER[p]')
TAG_PATH_NAME_MAP = {
    ('boundary', 'bulk_molecules_report', 'EG10040-MONOMER[p]'): 'AmpC',
    (
        'boundary', 'bulk_molecules_report', 'TRANS-CPLX-201[s]'
    ): 'AcrAB-TolC',
}
ENVIRONMENT_SECTION_FIELDS = ('GLC',)
ENVIRONMENT_SECTION_TIMES = (231, 4851, 9471, 13860, 18480, 23100)
AGENTS_TO_TRACE = (
    '0_wcecoli001100', '0_wcecoli0101011000', '0_wcecoli1001110001',
    '0_wcecoli1001110000', '0_wcecoli010101110011', '0_wcecoli011100',
    '0_wcecoli100001')
COLONY_MASS_PATH = ('mass',)
EXPRESSION_SURVIVAL_TIME_RANGE = (0.5, 1)
NUM_SNAPSHOTS = 5
FIG_OUT_DIR = os.path.join(OUT_DIR, 'figs')
FILE_EXTENSION = 'pdf'
EXPERIMENT_IDS = {
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
    'expression_survival': '20210301.030213',
    'death_snapshots': '20210301.030213',
    'phylogeny': '20210301.030213',
}
METADATA_FILE = 'metadata.json'
STATS_FILE = 'stats.json'


def exec_shell(tokens, timeout=10):
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


def get_metadata():
    '''Get information on which experiments and code were used.'''
    metadata = {
        'git_hash': exec_shell(['git', 'rev-parse', 'HEAD'])[0],
        'git_branch': exec_shell(
            ['git', 'symbolic-ref', '--short', 'HEAD'])[0],
        'git_status': exec_shell(
            ['git', 'status', '--porcelain'])[0].split('\n'),
        'time': datetime.utcnow().isoformat() + '+00:00',
        'python': sys.version.splitlines()[0],
        'experiment_ids': EXPERIMENT_IDS,
    }
    return metadata


def get_experiment_ids(id_obj):
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


def get_data(args, experiment_ids):
    '''Load all the data we'll need to generate the figures.'''
    unique_ids = set(experiment_ids)
    all_data = {}
    for experiment_id in unique_ids:
        all_data[experiment_id] = Analyzer.get_data(args, experiment_id)
    return all_data


def make_expression_heterogeneity_fig(
        data, environment_config, name_base):
    '''Figure shows heterogeneous expression within wcEcoli agents.'''
    tags_data = Analyzer.format_data_for_tags(data, environment_config)
    tagged_molecules = list(TAG_PATH_NAME_MAP.keys())
    colors = tuple(COLORS[name] for name in ('orange', 'cyan'))
    plot_config = {
        'out_dir': FIG_OUT_DIR,
        'tagged_molecules': tagged_molecules,
        'filename': '{}.{}'.format(name_base, FILE_EXTENSION),
        'tag_path_name_map': TAG_PATH_NAME_MAP,
        'tag_label_size': 48,
        'default_font_size': 48,
        'n_snapshots': NUM_SNAPSHOTS,
        'hues': {
            tag: mcolors.rgb_to_hsv(mcolors.to_rgb(colors[i]))[0]
            for i, tag in enumerate(tagged_molecules)
        },
    }
    plot_tags(tags_data, plot_config)


def make_expression_distributions_fig(replicates_raw_data):
    '''Figure shows the distributions of expression values.'''
    replicates_data = []
    colors = list(COLORS.values())
    for i, raw_data in enumerate(replicates_raw_data):
        color = colors[i]
        end_expression_table = raw_data_to_end_expression_table(
            raw_data,
            {val: key for key, val in TAG_PATH_NAME_MAP.items()})
        data = {
            key: end_expression_table[key]
            for key in end_expression_table.columns
            if key != VOLUME_KEY
        }
        replicates_data.append((data, color))
    fig = get_ridgeline_plot(
        replicates_data,
        point_alpha=1,
        overlap=-0.1,
        num_bins=100,
        jitter=0,
        x_label='Protein Concentration (counts/fL)',
        y_label='Distribution Density')
    fig.savefig(
        os.path.join(
            FIG_OUT_DIR,
            'expression_distributions.{}'.format(FILE_EXTENSION),
        )
    )


def make_snapshots_figure(
        data, environment_config, name, fields, agent_fill_color=None,
        agent_alpha=1):
    '''Make a figure of snapshots

    Parameters:
        data (dict): The experiment data.
        environment_config (dict): Environment parameters.
        name (str): Name of the output file (excluding file extension).
        fields (list): List of the names of fields to include.
    '''
    snapshots_data = Analyzer.format_data_for_snapshots(
        data, environment_config)
    if not fields:
        data = {
            key: val
            for key, val in data.items() if key != 'fields'
        }
    plot_config = {
        'out_dir': FIG_OUT_DIR,
        'filename': '{}.{}'.format(name, FILE_EXTENSION),
        'include_fields': fields,
        'field_label_size': 54,
        'default_font_size': 54,
        'agent_fill_color': agent_fill_color,
        'agent_alpha': agent_alpha,
        'n_snapshots': NUM_SNAPSHOTS,
    }
    plot_snapshots(snapshots_data, plot_config)


def make_growth_fig(basal_data, anaerobic_data):
    '''Make plot of colony mass of basal and anaerobic colonies.'''
    data_dict = {
        'basal': basal_data,
        'anaerobic': anaerobic_data,
    }
    fig, stats = get_total_mass_plot(data_dict, tuple(COLORS.values()))
    fig.savefig(os.path.join(
        FIG_OUT_DIR, 'growth.{}'.format(FILE_EXTENSION)))
    return stats


def make_threshold_scan_fig(data_dict):
    '''Plot colony mass curves with various antibiotic thresholds.'''
    fig, stats = get_total_mass_plot(data_dict, tuple(COLORS.values()))
    fig.savefig(os.path.join(
        FIG_OUT_DIR, 'threshold_scan.{}'.format(FILE_EXTENSION)))
    return stats


def make_expression_survival_fig(data, search_data):
    '''Make expression-survival figure.'''
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
        trace_agents=AGENTS_TO_TRACE,
    )
    fig.savefig(os.path.join(
        FIG_OUT_DIR, 'expression_survival.{}'.format(
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
    )
    fig.savefig(os.path.join(
        FIG_OUT_DIR, 'expression_survival_labeled.{}'.format(
            FILE_EXTENSION)
    ))


def make_expression_survival_dotplots(data):
    '''Make expression-survival dotplots.'''
    fig = plot_expression_survival_dotplot(
        data, PUMP_PATH, 'Final [AcrAB-TolC] (µM)',
        scaling=1e3,
        time_range=EXPRESSION_SURVIVAL_TIME_RANGE)
    fig.savefig(os.path.join(
        FIG_OUT_DIR, 'expression_survival_pump.{}'.format(
            FILE_EXTENSION)))
    fig = plot_expression_survival_dotplot(
        data, BETA_LACTAMASE_PATH, 'Final [AmpC] (µM)',
        scaling=1e3,
        time_range=EXPRESSION_SURVIVAL_TIME_RANGE)
    fig.savefig(os.path.join(
        FIG_OUT_DIR,
        'expression_survival_beta_lactamase.{}'.format(
            FILE_EXTENSION)))


def make_expression_survival_scan_fig(data, parameters, search_data):
    '''Plot expression-survival parameter scan figure.'''
    fig = plot_expression_survival_scan(
        data, parameters['agent_name'], '[AcrAB-TolC] (µM)',
        '[AmpC] (µM)', search_data['x_values'], search_data['y_values'],
        search_data['precision'], scaling=1e3)
    fig.savefig(os.path.join(
        FIG_OUT_DIR,
        'expression_survival_scan.{}'.format(FILE_EXTENSION)
    ))


def make_environment_section(data, base_name):
    '''Plot field concentrations in cross-section of final enviro.'''
    t_final = max(data[0].keys())
    fields_ts = []
    section_times = [
        float(time) for time in ENVIRONMENT_SECTION_TIMES]
    for i, replicate in enumerate(data):
        fields_ts.append(dict())
        for time in section_times:
            fields_ts[i][time] = {
                name: field
                for name, field in get_in(
                    replicate[time], FIELDS_PATH).items()
                if name in ENVIRONMENT_SECTION_FIELDS
            }
    bounds = get_in(data[0][t_final], BOUNDS_PATH)
    fig = get_enviro_sections_plot(fields_ts, bounds,
        section_location=0.5)
    fig.savefig(
        os.path.join(FIG_OUT_DIR, '{}.{}'.format(
            base_name, FILE_EXTENSION)))


def make_phylogeny_plot(data):
    '''Plot phylogenetic tree'''
    plot_phylogeny(data, os.path.join(
        FIG_OUT_DIR, 'phylogeny.{}').format(FILE_EXTENSION),
        time_range=EXPRESSION_SURVIVAL_TIME_RANGE)


def main():
    '''Generate all figures.'''
    if not os.path.exists(FIG_OUT_DIR):
        os.makedirs(FIG_OUT_DIR)
    with open(os.path.join(FIG_OUT_DIR, METADATA_FILE), 'w') as f:
        json.dump(get_metadata(), f, indent=4)
    stats = {}
    parser = argparse.ArgumentParser()
    Analyzer.add_connection_args(parser)
    parser.add_argument(
        'scan_data', type=str, help='Path to parameter scan data.')
    parser.add_argument(
        'search_data', type=str, help='Path to boundary search data.')
    args = parser.parse_args()

    experiment_ids = get_experiment_ids(EXPERIMENT_IDS)
    all_data = get_data(args, experiment_ids)

    expression_distribution_data = []
    for experiment_id in EXPERIMENT_IDS['expression_distributions']:
        data, _ = all_data[experiment_id]
        expression_distribution_data.append(data)
    make_expression_distributions_fig(expression_distribution_data)

    for i, experiment_id in enumerate(
            EXPERIMENT_IDS['expression_heterogeneity']):
        make_expression_heterogeneity_fig(
            *all_data[experiment_id],
            'expression_heterogeneity_{}'.format(i))

    for i, experiment_id in enumerate(EXPERIMENT_IDS['growth_basal']):
        make_snapshots_figure(
            *all_data[experiment_id], 'growth_basal_{}'.format(i), [])

    for i, experiment_id in enumerate(
            EXPERIMENT_IDS['growth_anaerobic']):
        make_snapshots_figure(
            *all_data[experiment_id], 'growth_anaerobic_{}'.format(i),
            [])

    basal_data = [
        all_data[experiment_id][0]
        for experiment_id in EXPERIMENT_IDS['growth_basal']
    ]
    anaerobic_data = [
        all_data[experiment_id][0]
        for experiment_id in EXPERIMENT_IDS['growth_anaerobic']
    ]
    stats['growth_fig'] = make_growth_fig(basal_data, anaerobic_data)

    for i, experiment_id in enumerate(
            EXPERIMENT_IDS['enviro_heterogeneity']):
        make_snapshots_figure(
            *all_data[experiment_id],
            'enviro_heterogeneity_{}'.format(i), ['GLC'], 'black', 0.2)

    enviro_section_data = []
    for i, experiment_id in enumerate(
            EXPERIMENT_IDS['enviro_section']):
        enviro_section_data.append(all_data[experiment_id][0])
    make_environment_section(
        enviro_section_data, 'enviro_section')

    data_dict = dict()
    for key, exp_ids in EXPERIMENT_IDS['threshold_scan'].items():
        data_dict[key] = [all_data[exp_id][0] for exp_id in exp_ids]
    stats['threshold_scan'] = make_threshold_scan_fig(data_dict)

    with open(args.search_data, 'r') as f:
        search_data = json.load(f)

    make_expression_survival_fig(
        all_data[EXPERIMENT_IDS['expression_survival']][0], search_data)

    make_expression_survival_dotplots(
        all_data[EXPERIMENT_IDS['expression_survival']][0])

    make_expression_survival_scan_fig(
        *load_scan_data(args.scan_data), search_data)

    make_phylogeny_plot(
        all_data[EXPERIMENT_IDS['phylogeny']][0])

    make_snapshots_figure(
        *all_data[EXPERIMENT_IDS['death_snapshots']],
        'death_snapshots', ['nitrocefin'], 'green')

    with open(os.path.join(FIG_OUT_DIR, STATS_FILE), 'w') as f:
        json.dump(stats, f, indent=4)


if __name__ == '__main__':
    main()
