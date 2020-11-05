'''Generate development analysis plots for colony simulations

For usage information, execute:

    python analyze.py -h
'''

import argparse
import os

from vivarium_cell.analysis.analyze import Analyzer
from vivarium.core.experiment import get_in

from src.constants import FIELDS_PATH, BOUNDS_PATH
from src.environment_cross_sections import get_enviro_sections_plot
from src.phylogeny import plot_phylogeny


#: Configuration dictionary for multigen agent timeseries plot.
TIMESERIES_CONFIG = {
    'skip_paths': [
        ('boundary', 'wcecoli_fields_null'),
    ],
}
#: Configuration dictionary for snapshots plot.
SNAPSHOTS_CONFIG = {
    'include_fields': ['nitrocefin', 'GLC'],
    'field_label_size': 54,
    'default_font_size': 54,
}
#: Configuration dictionary for tagged molecules plot.
TAGS_CONFIG = {
    'tag_label_size': 54,
    'default_font_size': 54,
}
#: Fields to Plot in Environment Cross-Section
ENVIRONMENT_SECTION_FIELDS = None


class ColonyAnalyzer(Analyzer):
    '''Analyze a wcEcoli colony simulation'''

    def _get_parser(self) -> argparse.ArgumentParser:
        '''
        Let the parser handle requests to generate phylogeny tree.
        '''
        parser = super()._get_parser()
        parser.add_argument(
            '--phylogeny', '-y',
            action='store_true',
            default=False,
            help='Plot agent phylogeny',
        )
        parser.add_argument(
            '--environment_section', '-e',
            type=str,
            help=(
                'Plot cross-sections of environmental fields. '
                'Specify "flat" or "mid".'
            ),
        )
        return parser

    def plot(self, args: argparse.Namespace) -> None:
        '''Let analyzer generate phylogeny plots.'''
        super().plot(args)
        if args.phylogeny:
            plot_phylogeny(
                self.data,
                os.path.join(self.out_dir, 'phylogeny.png'),
            )
        if args.environment_section:
            assert args.environment_section in ('flat', 'mid')
            t_final = max(self.data.keys())
            fields = get_in(self.data[t_final], FIELDS_PATH)
            if ENVIRONMENT_SECTION_FIELDS:
                fields = {
                    key: val for key, val in fields.items()
                    if key in ENVIRONMENT_SECTION_FIELDS
                }
            bounds = get_in(self.data[t_final], BOUNDS_PATH)
            flat_bins = args.environment_section == 'flat'
            fig = get_enviro_sections_plot(fields, bounds,
                    section_location=0.4, flat_bins=flat_bins)
            fig.savefig(
                os.path.join(self.out_dir, 'enviro_sections.png'))


def main() -> None:
    '''Handle CLI args and generate plots.'''
    analyzer = ColonyAnalyzer(
        timeseries_config=TIMESERIES_CONFIG,
        snapshots_config=SNAPSHOTS_CONFIG,
    )
    analyzer.run()


if __name__ == '__main__':
    main()
