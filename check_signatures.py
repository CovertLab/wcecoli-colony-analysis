'''Check that all commits are signed by an allowed committer.'''


from typing import Sequence, Tuple

import os
import subprocess
import sys


TIMEOUT = 60 * 10
ALLOWED_COMMITTERS_PATH = os.path.join(
    os.path.dirname(__file__), 'ALLOWED_COMMITTERS')


def exec_shell(
        tokens: Sequence[str],
        timeout: int = 10,
        ) -> Tuple[str, str]:
    '''Execute a shell command and return the output.

    Args:
        tokens: Sequence of space-separated tokens that will be passed
            to the shell.
        timeout: Seconds to wait for the command to complete before
            abording.

    Returns:
        Tuple of the stdout and stderr from the shell command.
    '''
    proc = subprocess.run(
        tokens,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
        env=None,
        universal_newlines=True,
        timeout=timeout)
    return proc.stdout.rstrip(), proc.stderr.rstrip()


def main() -> bool:
    '''Check that all commits were made by allowed committers.

    Returns:
        True if all commits are correctly signed, else False.
    '''
    # The ALLOWED_COMMITTERS file lists the fingerprints of keys allowed
    # to sign commits, one per line.
    with open(ALLOWED_COMMITTERS_PATH, 'r') as f:
        allowed_committers = {line.rstrip() for line in f}
    # Run `git log` with a format such that each commit is shown on its
    # own line as the commit hash, a space, and the fingerprint of the
    # primary key that signed the commit.
    log, _ = exec_shell(['git', 'log', '--format=%H %GP'], TIMEOUT)
    log_lines = log.split('\n')
    bad_commits = []  # List of commits not correctly signed
    for line in log_lines:
        split = line.split()  # Split a line by the space delimiter
        if len(split) == 1:
            # Only the commit was found, so the commit was not signed.
            bad_commits.append(split[0])
            continue
        assert len(split) == 2
        commit, committer = split
        # Check that the committer fingerprint is in the set of allowed
        # committers.
        if committer not in allowed_committers:
            bad_commits.append(commit)
    if bad_commits:
        print('Commits not signed by an allowed committer:')
        # Print the commits that were not correctly signed
        for commit in bad_commits:
            print('    ' + commit)
        return False
    return True


if __name__ == '__main__':
    # Exit with code 0 if main() returns True, else exit with code 1.
    sys.exit(0 if main() else 1)
