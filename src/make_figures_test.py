import os
from unittest.mock import patch

from src import make_figures


EXPERIMENT_IDS = {
    'expression_distributions': (
        '20210125.180927', '20210125.182741', '20210125.184150'),
    'expression_heterogeneity': (
        '20210125.180927', '20210125.182741', '20210125.184150'),
    'enviro_heterogeneity': (
        '20210125.180927', '20210125.182741', '20210125.184150'),
    'enviro_section': (
        '20210125.180927', '20210125.182741', '20210125.184150'),
    'growth_basal': (
        '20210125.180927', '20210125.182741', '20210125.184150'),
    'growth_anaerobic': (
        '20210125.181216', '20210125.184703', '20210125.185045'),
    'threshold_scan': {
        '0.01 mM': ('20210125.182021',),
        '0.02 mM': ('20210125.182506',),
        '0.03 mM': ('20210125.183030',),
        '0.04 mM': ('20210125.183244',),
    },
    'expression_survival': '20210125.182506',
    'death_snapshots': '20210125.182506',
    'phylogeny': '20210125.182506',
}
ENVIRONMENT_SECTION_TIMES = (1, 21, 41, 60, 80, 100)
AGENTS_TO_TRACE = ('0_wcecoli',)
FIG_OUT_DIR = os.path.join(make_figures.OUT_DIR, 'figs_test')


def main():
    make_figures.EXPERIMENT_IDS = EXPERIMENT_IDS
    make_figures.FIG_OUT_DIR = FIG_OUT_DIR
    make_figures.ENVIRONMENT_SECTION_TIMES = ENVIRONMENT_SECTION_TIMES
    make_figures.AGENTS_TO_TRACE = AGENTS_TO_TRACE
    make_figures.main()


if __name__ == '__main__':
    main()
