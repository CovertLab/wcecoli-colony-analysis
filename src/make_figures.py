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
from vivarium_cell.plots.expression_survival_dotplot import (
    plot_expression_survival)

from src.constants import OUT_DIR, FIELDS_PATH, BOUNDS_PATH
from src.total_mass import get_total_mass_plot
from src.environment_cross_sections import get_enviro_sections_plot


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
ENVIRONMENT_SECTION_FIELDS = (
    'GLC', 'AMMONIUM', 'PI', 'HYPOXANTHINE', 'K+', 'TRP',
    'ASN', 'L_ALPHA_ALANINE', 'SULFATE', 'PROTON')
COLONY_MASS_PATH = ('mass',)
FIG_OUT_DIR = os.path.join(OUT_DIR, 'figs')
FILE_EXTENSION = 'pdf'
EXPRESSION_HETEROGENEITY_ID = '20200820.202016'
ENVIRO_HETEROGENEITY_ID = '20200824.165625'
ENVIRO_SECTION_ID = ENVIRO_HETEROGENEITY_ID
GROWTH_BASAL_ID = EXPRESSION_HETEROGENEITY_ID
GROWTH_ANAEROBIC_ID = '20200820.235622'
THRESHOLD_SCAN_IDS = {
    '0.002 mM': '20200818.174841',
    '0.01775 mM': '20200819.175108',
    '0.018875 mM': '20200819.203802',
    '0.02 mM': '20200817.224609',
    '0.025 mM': '20200821.172829',
    '0.026 mM': '20200824.141256',
    '0.0275 mM': '20200823.195457',
    '0.03 mM': '20200821.142922',
}
EXPRESSION_SURVIVAL_ID = '20200817.224609'
PUMP_TIMESERIES_ID = EXPRESSION_SURVIVAL_ID
EXPRESSION_SURVIVAL_TIME_RANGE = (0.5, 1)
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
        'pump_timeseries_id': PUMP_TIMESERIES_ID,
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
        data, PUMP_PATH,
        'Average AcrAB-TolC Concentration (mmol/L) Over Cell Lifetime',
        EXPRESSION_SURVIVAL_TIME_RANGE,
    )
    fig.savefig(os.path.join(
        FIG_OUT_DIR, 'expression_survival_acrABtolC.{}'.format(
            FILE_EXTENSION)
    ))

    fig = plot_expression_survival(
        data, BETA_LACTAMASE_PATH,
        'Average AmpC Concentration (mmol/L) Over Cell Lifetime',
        EXPRESSION_SURVIVAL_TIME_RANGE,
    )
    fig.savefig(os.path.join(
        FIG_OUT_DIR, 'expression_survival_ampC.{}'.format(
            FILE_EXTENSION)
    ))


def make_pump_timeseries_fig(data):
    '''Plot concentrations over time for all agents.'''
    settings = {
        'include_paths': [PUMP_PATH],
        'titles_map': {
            PUMP_PATH: 'AcrAB-TolC Concentration',
        },
        'ylabels_map': {
            PUMP_PATH: 'Concentration (mM)',
        },
    }
    plot_agents_multigen(
        data, settings, FIG_OUT_DIR,
        'acrABtolC_timeseries.{}'.format(FILE_EXTENSION)
    )
    settings = {
        'include_paths': [BETA_LACTAMASE_PATH],
        'titles_map': {
            BETA_LACTAMASE_PATH: 'AmpC Concentration',
        },
        'ylabels_map': {
            BETA_LACTAMASE_PATH: 'Concentration (mM)',
        },
    }
    plot_agents_multigen(
        data, settings, FIG_OUT_DIR,
        'ampC_timeseries.{}'.format(FILE_EXTENSION)
    )


def make_environment_section(data):
    '''Plot field concentrations in cross-section of final enviro.'''
    t_final = max(data.keys())
    fields = get_in(data[t_final], FIELDS_PATH)
    fields = {
        key: val for key, val in fields.items()
        if key in ENVIRONMENT_SECTION_FIELDS
    }
    bounds = get_in(data[t_final], BOUNDS_PATH)
    fig = get_enviro_sections_plot(fields, bounds,
            section_location=0.5, flat_bins=False)
    fig.savefig(
        os.path.join(FIG_OUT_DIR, 'enviro_sections.{}'.format(
            FILE_EXTENSION))
    )


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
        data, environment_config, 'growth_basal', ['nitrocefin'])

    data, environment_config = Analyzer.get_data(
        args, GROWTH_ANAEROBIC_ID)
    make_snapshots_figure(
        data, environment_config, 'growth_anaerobic', ['nitrocefin'])

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

    data, _ = Analyzer.get_data(args, EXPRESSION_SURVIVAL_ID)
    make_expression_survival_fig(data)

    if EXPRESSION_SURVIVAL_ID != PUMP_TIMESERIES_ID:
        data, _ = Analyzer.get_data(args, PUMP_TIMESERIES_ID)
    make_pump_timeseries_fig(data)


if __name__ == '__main__':
    main()
