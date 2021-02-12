'''Archive all the simulations used to generate figures.'''

import os
import shutil
import sys
import tempfile

from tqdm import tqdm

from src.make_figures import (
    EXPERIMENT_IDS, get_experiment_ids, exec_shell)


OUT_PATH = 'archived_simulations.tar.gz'
TIMEOUT = 10 * 60  # seconds
ARCHIVES_FOLDER = 'archived_simulations'


def main():
    experiment_ids = set(get_experiment_ids(EXPERIMENT_IDS))

    with tempfile.TemporaryDirectory() as archive_dir:
        os.mkdir(os.path.join(archive_dir, ARCHIVES_FOLDER))
        print('Downloading Experiments')
        for experiment_id in tqdm(experiment_ids):
            archive_file = '{}.json'.format(experiment_id)
            if not os.path.exists(archive_file):
                args = sys.argv[1:] + ['download', experiment_id]
                exec_shell(
                    ['python', '-m', 'scripts.access_db'] + args,
                    timeout=TIMEOUT)
        print('Moving archives to temporary directory', archive_dir)
        for experiment_id in experiment_ids:
            archive_file = '{}.json'.format(experiment_id)
            shutil.move(
                archive_file,
                os.path.join(archive_dir, ARCHIVES_FOLDER, archive_file)
            )
        print('Creating compressed archive', OUT_PATH)
        exec_shell([
            'tar', '-C', archive_dir, '-c', '-z', '-f', OUT_PATH,
            ARCHIVES_FOLDER])


if __name__ == '__main__':
    main()
