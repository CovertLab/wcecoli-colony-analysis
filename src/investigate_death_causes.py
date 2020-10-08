'''Script for investigating what characteristics predict death.
'''


import os

from vivarium.plots.agents_multigen import plot_agents_multigen
from vivarium_cell.analysis.analyze import Analyzer
from vivarium_cell.plots.expression_survival_dotplot import (
    plot_expression_survival)
from vivarium_cell.plots.multibody_physics import plot_tags

from src.constants import (
    ACRAB_TOLC_KEY,
    BETA_LACTAMASE_KEY,
)
from src.investigate_utils import (
    parse_args_retrieve_data_by_id,
    filter_raw_data_by_time,
    split_raw_data_by_survival,
)


PUMP_PATH = (
    'boundary', 'bulk_molecule_concentrations', ACRAB_TOLC_KEY)
BETA_LACTAMASE_PATH = (
    'boundary', 'bulk_molecule_concentrations', BETA_LACTAMASE_KEY)
ANTIBIOTIC_TIME_RANGE = (0.5, 1)


def main() -> None:
    '''Run analyses'''
    data, environment_config, out_dir = parse_args_retrieve_data_by_id()

    fig_pump = plot_expression_survival(
        data, PUMP_PATH,
        'Average [AcrAB-TolC] (mmol/L) over Cell Lifetime',
        ANTIBIOTIC_TIME_RANGE,
    )
    fig_pump.savefig(os.path.join(out_dir, 'pump'))
    fig_beta_lactamase = plot_expression_survival(
        data, BETA_LACTAMASE_PATH,
        'Average [Beta-Lactamase] (mmol/L) over Cell Lifetime',
        ANTIBIOTIC_TIME_RANGE,
    )
    fig_beta_lactamase.savefig(os.path.join(out_dir, 'beta_lactamase'))

    multigen_settings = {
        'include_paths': [
            PUMP_PATH,
            BETA_LACTAMASE_PATH,
            ('boundary', 'bulk_molecules_report', ACRAB_TOLC_KEY),
            ('boundary', 'bulk_molecules_report', BETA_LACTAMASE_KEY),
            ('boundary', 'cytoplasm', 'nitrocefin_hydrolyzed'),
            ('boundary', 'cytoplasm', 'nitrocefin'),
            ('boundary', 'external', 'nitrocefin'),
            ('boundary', 'dead'),
        ],
    }
    filtered = filter_raw_data_by_time(data, ANTIBIOTIC_TIME_RANGE)
    survive_data, die_data = split_raw_data_by_survival(filtered)
    plot_agents_multigen(
        survive_data, multigen_settings, out_dir, 'survive')
    plot_agents_multigen(
        die_data, multigen_settings, out_dir, 'die')

    tag_path_name_map = {
        PUMP_PATH: 'AcrAB-TOlC',
        BETA_LACTAMASE_PATH: 'Beta-Lactamase',
    }
    tags_config = {
        'out_dir': out_dir,
        'tagged_molecules': tag_path_name_map.keys(),
        'filename': 'expression',
        'tag_path_name_map': tag_path_name_map,
    }
    tags_data = Analyzer.format_data_for_tags(data, environment_config)
    plot_tags(tags_data, tags_config)


if __name__ == '__main__':
    main()
