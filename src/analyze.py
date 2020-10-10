'''Generate development analysis plots for colony simulations

For usage information, execute:

    python analyze.py -h
'''

import argparse
import os

from vivarium_cell.analysis.analyze import Analyzer

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
        return parser

    def plot(self, args: argparse.Namespace) -> None:
        '''Let analyzer generate phylogeny plots.'''
        super().plot(args)
        if args.phylogeny:
            plot_phylogeny(
                self.data,
                os.path.join(self.out_dir, 'phylogeny.png'),
            )


def main() -> None:
    '''Handle CLI args and generate plots.'''
    analyzer = ColonyAnalyzer(
        timeseries_config=TIMESERIES_CONFIG,
        snapshots_config=SNAPSHOTS_CONFIG,
    )
    analyzer.run()


if __name__ == '__main__':
    main()
