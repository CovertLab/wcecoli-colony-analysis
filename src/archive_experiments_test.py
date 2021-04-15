from src import archive_experiments
from src.make_figures_test import EXPERIMENT_IDS


def main() -> None:
    '''Archive experiments after patching IDs constant.'''
    archive_experiments.EXPERIMENT_IDS = EXPERIMENT_IDS
    archive_experiments.main()


if __name__ == '__main__':
    main()
