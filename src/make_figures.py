'''Generate figures for wcEcoli colony simulation

For usage information, run:
    python make_figures.py -h
'''
import argparse
import json
import os
import sys

from vivarium.plots.agents_multigen import plot_agents_multigen
from vivarium.core.experiment import get_in
from vivarium_cell.analysis.analyze import Analyzer
from vivarium_cell.plots.multibody_physics import (
    plot_snapshots,
    plot_tags,
)
from src.expression_survival import plot_expression_survival

from src.constants import OUT_DIR, FIELDS_PATH, BOUNDS_PATH
from src.total_mass import get_total_mass_plot
from src.environment_cross_sections import get_enviro_sections_plot
from src.phylogeny import plot_phylogeny


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
COLONY_MASS_PATH = ('mass',)
FIG_OUT_DIR = os.path.join(OUT_DIR, 'figs')
FILE_EXTENSION = 'pdf'
EXPRESSION_HETEROGENEITY_ID = '20201119.150828'
ENVIRO_HETEROGENEITY_ID = '20201119.150828'
ENVIRO_SECTION_ID = ENVIRO_HETEROGENEITY_ID
GROWTH_BASAL_ID = EXPRESSION_HETEROGENEITY_ID
GROWTH_ANAEROBIC_ID = '20201221.194828'
THRESHOLD_SCAN_IDS = {
    '0.01 mM': '20201228.172246',
    '0.02 mM': '20201228.211700',
    '0.03 mM': '20201229.160649',
    '0.04 mM': '20201230.191552',
}
EXPRESSION_SURVIVAL_ID = THRESHOLD_SCAN_IDS['0.02 mM']
EXPRESSION_SURVIVAL_TIME_RANGE = (0.5, 1)
DEATH_SNAPSHOTS_ID = THRESHOLD_SCAN_IDS['0.02 mM']
PHYLOGENY_ID = THRESHOLD_SCAN_IDS['0.02 mM']
METADATA_FILE = 'metadata.json'


def get_metadata():
    '''Get information on which experiments and code were used.'''
    metadata = {
        #'git_hash': fp.run_cmdline('git rev-parse HEAD'),
        #'git_branch': fp.run_cmdline('git symbolic-ref --short HEAD'),
        #'time': fp.timestamp(),
        'python': sys.version.splitlines()[0],
        'expression_heterogeneity_id': EXPRESSION_HETEROGENEITY_ID,
        'enviro_heterogeneity_id': ENVIRO_HETEROGENEITY_ID,
        'enviro_section_id': ENVIRO_SECTION_ID,
        'growth_basal_id': GROWTH_BASAL_ID,
        'growth_anaerobic_id': GROWTH_ANAEROBIC_ID,
        'threshold_scan_ids': THRESHOLD_SCAN_IDS,
        'expression_survival_id': EXPRESSION_SURVIVAL_ID,
        'death_snapshots_id': DEATH_SNAPSHOTS_ID,
        'phylogeny_id': PHYLOGENY_ID,
    }
    return metadata


def make_expression_heterogeneity_fig(data, environment_config):
    '''Figure shows heterogeneous expression within wcEcoli agents.'''
    tags_data = Analyzer.format_data_for_tags(data, environment_config)
    plot_config = {
        'out_dir': FIG_OUT_DIR,
        'tagged_molecules': TAG_PATH_NAME_MAP.keys(),
        'filename': 'expression_heterogeneity.{}'.format(
            FILE_EXTENSION),
        'tag_path_name_map': TAG_PATH_NAME_MAP,
        'tag_label_size': 48,
        'default_font_size': 48,
    }
    plot_tags(tags_data, plot_config)


def make_snapshots_figure(data, environment_config, name, fields):
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
        data = {key: val for key, val in data.items if key != 'fields'}
    plot_config = {
        'out_dir': FIG_OUT_DIR,
        'filename': '{}.{}'.format(name, FILE_EXTENSION),
        'include_fields': fields,
        'field_label_size': 54,
        'default_font_size': 54,
    }
    plot_snapshots(snapshots_data, plot_config)


def make_growth_fig(basal_data, anaerobic_data):
    '''Make plot of colony mass of basal and anaerobic colonies.'''
    data_dict = {
        'basal': basal_data,
        'anaerobic': anaerobic_data,
    }
    fig = get_total_mass_plot(data_dict)
    fig.savefig(os.path.join(
        FIG_OUT_DIR, 'growth.{}'.format(FILE_EXTENSION)))


def make_threshold_scan_fig(data_dict):
    '''Plot colony mass curves with various antibiotic thresholds.'''
    fig = get_total_mass_plot(data_dict)
    fig.savefig(os.path.join(
        FIG_OUT_DIR, 'threshold_scan.{}'.format(FILE_EXTENSION)))


def make_expression_survival_fig(data):
    '''Make expression-survival dotplot figure.'''
    fig = plot_expression_survival(
        data, PUMP_PATH, BETA_LACTAMASE_PATH,
        'Average AcrAB-TolC Concentration (mmol/L) Over Cell Lifetime',
        'Average AmpC Concentration (mmol/L) Over Cell Lifetime',
        EXPRESSION_SURVIVAL_TIME_RANGE,
    )
    fig.savefig(os.path.join(
        FIG_OUT_DIR, 'expression_survival.{}'.format(
            FILE_EXTENSION)
    ))


def make_environment_section(data):
    '''Plot field concentrations in cross-section of final enviro.'''
    t_final = max(data.keys())
    fields_ts = dict()
    section_times = [
        float(time) for time in ENVIRONMENT_SECTION_TIMES]
    for time in section_times:
        fields_ts[time] = {
            name: field
            for name, field in get_in(
                data[time], FIELDS_PATH).items()
            if name in ENVIRONMENT_SECTION_FIELDS
        }
    bounds = get_in(data[t_final], BOUNDS_PATH)
    fig = get_enviro_sections_plot(fields_ts, bounds,
            section_location=0.5, flat_bins=False)
    fig.savefig(
        os.path.join(FIG_OUT_DIR, 'enviro_sections.{}'.format(
            FILE_EXTENSION)))


def make_phylogeny_plot(data):
    '''Plot phylogenetic tree'''
    plot_phylogeny(data, os.path.join(
        FIG_OUT_DIR, 'phylogeny.{}').format(FILE_EXTENSION))


def main():
    '''Generate all figures.'''
    if not os.path.exists(FIG_OUT_DIR):
        os.makedirs(FIG_OUT_DIR)
    with open(os.path.join(FIG_OUT_DIR, METADATA_FILE), 'w') as f:
        json.dump(get_metadata(), f)
    parser = argparse.ArgumentParser()
    Analyzer.add_connection_args(parser)
    args = parser.parse_args()

    data, environment_config = Analyzer.get_data(
        args, EXPRESSION_HETEROGENEITY_ID)
    make_expression_heterogeneity_fig(data, environment_config)

    if GROWTH_BASAL_ID != EXPRESSION_HETEROGENEITY_ID:
        data, environment_config = Analyzer.get_data(
            args, GROWTH_BASAL_ID)
    data_growth_basal = data
    make_snapshots_figure(
        data, environment_config, 'growth_basal', [])

    data, environment_config = Analyzer.get_data(
        args, GROWTH_ANAEROBIC_ID)
    make_snapshots_figure(
        data, environment_config, 'growth_anaerobic', [])

    make_growth_fig(data_growth_basal, data)
    del data_growth_basal

    if ENVIRO_HETEROGENEITY_ID != GROWTH_ANAEROBIC_ID:
        data, environment_config = Analyzer.get_data(
            args, ENVIRO_HETEROGENEITY_ID)
    make_snapshots_figure(
        data, environment_config, 'enviro_heterogeneity', ['GLC'])

    if ENVIRO_SECTION_ID != ENVIRO_HETEROGENEITY_ID:
        data, environment_config = Analyzer.get_data(
            args, ENVIRO_SECTION_ID)
    make_environment_section(data)

    data_dict = dict()
    for key, exp_id in THRESHOLD_SCAN_IDS.items():
        exp_data, _ = Analyzer.get_data(args, exp_id)
        data_dict[key] = exp_data
    make_threshold_scan_fig(data_dict)
    del data_dict

    data, environment_config = Analyzer.get_data(
        args, EXPRESSION_SURVIVAL_ID)
    make_expression_survival_fig(data)

    if PHYLOGENY_ID != EXPRESSION_SURVIVAL_ID:
        data, environment_config = Analyzer.get_data(args, PHYLOGENY_ID)
    make_phylogeny_plot(data)

    if DEATH_SNAPSHOTS_ID != PHYLOGENY_ID:
        data, environment_config = Analyzer.get_data(
            args, DEATH_SNAPSHOTS_ID)
    make_snapshots_figure(
        data, environment_config, 'death_snapshots', ['nitrocefin'])


if __name__ == '__main__':
    main()
