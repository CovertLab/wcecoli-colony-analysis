from matplotlib import pyplot as plt
import numpy as np
from vivarium.library.topology import get_in

from src.constants import AGENTS_PATH


MONOMER_COUNTS_PATH = ('listeners', 'monomer_counts')
SUBMASS_PATHS = {
    'Total Dry Mass': ('listeners', 'mass', 'dry_mass'),
    'Protein Mass': ('listeners', 'mass', 'proteinMass'),
    'tRNA Mass': ('listeners', 'mass', 'tRnaMass'),
    'rRNA Mass': ('listeners', 'mass', 'rRnaMass'),
    'mRNA Mass': ('listeners', 'mass', 'mRnaMass'),
    'DNA Mass': ('listeners', 'mass', 'dnaMass'),
    'Small Molecule Mass': ('listeners', 'mass', 'smallMoleculeMass'),
}


def _get_timeseries(datasets, variable_path, agent_target):
    times = []
    values = []
    for dataset in datasets:
        dataset_times = []
        dataset_values = []
        for time, timepoint_data in dataset.items():
            if agent_target:
                agent_data = get_in(
                    timepoint_data, AGENTS_PATH + (agent_target,))
                if not agent_data:
                    # The target agent is not present.
                    continue
            else:
                agent_data = timepoint_data
            dataset_times.append(time)
            dataset_values.append(get_in(agent_data, variable_path))
        times.append(np.array(dataset_times))
        values.append(np.array(dataset_values))
    # Each value in values is an array with times on axis 0. If each
    # variable value is an array, then that array's indexes are on axis
    # 1.
    return times, values


def plot_proteome_comparison(ax, data, vivarium_agent='', wcecoli_agent=''):
    vivarium_data = data['vivarium-ecoli']
    wcecoli_data = data['wcecoli']

    _, vivarium_timeseries = _get_timeseries(
        vivarium_data, MONOMER_COUNTS_PATH, vivarium_agent)
    vivarium_avg_protein = np.array([
        timeseries.mean(axis=0)
        for timeseries in vivarium_timeseries
    ])
    _, wcecoli_timeseries = _get_timeseries(
        wcecoli_data, MONOMER_COUNTS_PATH, wcecoli_agent)
    wcecoli_avg_protein = np.array([
        timeseries.mean(axis=0)
        for timeseries in wcecoli_timeseries
    ])

    vivarium_avg_protein = np.log10(vivarium_avg_protein + 1)
    wcecoli_avg_protein = np.log10(wcecoli_avg_protein + 1)

    ax.set_xlabel('log10(wcEcoli protein counts + 1)')
    ax.set_ylabel('log10(vivarium-ecoli protein counts + 1)')

    (
        wcecoli_avg_protein_q25,
        wcecoli_avg_protein_q50,
        wcecoli_avg_protein_q75,
    ) = np.percentile(wcecoli_avg_protein, [25, 50, 75], axis=0)
    (
        vivarium_avg_protein_q25,
        vivarium_avg_protein_q50,
        vivarium_avg_protein_q75,
    ) = np.percentile(vivarium_avg_protein, [25, 50, 75], axis=0)
    ax.errorbar(
        wcecoli_avg_protein_q50,
        vivarium_avg_protein_q50,
        xerr=[
            wcecoli_avg_protein_q50 - wcecoli_avg_protein_q25,
            wcecoli_avg_protein_q75 - wcecoli_avg_protein_q50,
        ],
        yerr=[
            vivarium_avg_protein_q50 - vivarium_avg_protein_q25,
            vivarium_avg_protein_q75 - vivarium_avg_protein_q50,
        ],
        fmt='none',
        ecolor='black',
        elinewidth=1,
        alpha=0.2,
    )

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    lower = max(xmin, ymin)
    upper = min(xmax, ymax)
    ax.plot(
        [lower, upper], [lower, upper], linestyle='--', color='gray',
        zorder=-1, linewidth=1,
    )


def get_proteome_comparison_plot(
        data, vivarium_agent='', wcecoli_agent=''):
    fig, ax = plt.subplots()
    plot_proteome_comparison(ax, data, vivarium_agent, wcecoli_agent)
    return fig


def plot_submass_comparison(
        ax, data, colors, vivarium_agent='', wcecoli_agent=''):
    vivarium_data = data['vivarium-ecoli']
    wcecoli_data = data['wcecoli']

    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Submasses Normalized by Initial Masses')

    for color, (label, path) in zip(colors, SUBMASS_PATHS.items()):
        vivarium_times, vivarium_timeseries = _get_timeseries(
            vivarium_data, path, vivarium_agent)
        wcecoli_times, wcecoli_timeseries = _get_timeseries(
            wcecoli_data, path, wcecoli_agent)

        # Make sure all vivarium datasets are the same length.
        n_vivarium_timepoints = min([
            len(times) for times in vivarium_times])
        vivarium_time_arrays = [
            times[:n_vivarium_timepoints]
            for times in vivarium_times
        ]
        vivarium_timeseries = np.array([
            timeseries[:n_vivarium_timepoints]
            for timeseries in vivarium_timeseries
        ])

        # Make sure all wcecoli datasets are the same length.
        n_wcecoli_timepoints = min([
            len(times) for times in wcecoli_times])
        wcecoli_time_arrays = [
            times[:n_wcecoli_timepoints]
            for times in wcecoli_times
        ]
        wcecoli_timeseries = np.array([
            timeseries[:n_wcecoli_timepoints]
            for timeseries in wcecoli_timeseries
        ])

        vivarium_times = vivarium_time_arrays[0]
        for time_array in vivarium_time_arrays:
            assert np.all(vivarium_times == time_array)
        wcecoli_times = wcecoli_time_arrays[0]
        for time_array in wcecoli_time_arrays:
            assert np.all(wcecoli_times == time_array)

        vivarium_timeseries /= vivarium_timeseries[:, 0].reshape(
            (vivarium_timeseries.shape[0], 1))
        wcecoli_timeseries /= wcecoli_timeseries[:, 0].reshape(
            (wcecoli_timeseries.shape[0], 1))

        (
            wcecoli_q25,
            wcecoli_q50,
            wcecoli_q75,
        ) = np.percentile(wcecoli_timeseries, [25, 50, 75], axis=0)
        (
            vivarium_q25,
            vivarium_q50,
            vivarium_q75,
        ) = np.percentile(vivarium_timeseries, [25, 50, 75], axis=0)

        ax.fill_between(
            vivarium_times, vivarium_q25, vivarium_q75, color=color,
            alpha=0.2)
        ax.fill_between(
            wcecoli_times, wcecoli_q25, wcecoli_q75, color=color,
            alpha=0.2)
        ax.plot(
            vivarium_times, vivarium_q50,
            label=f'{label} (vivarium-ecoli)', linestyle='-',
            linewidth=1, color=color)
        ax.plot(
            wcecoli_times, wcecoli_q50, label=f'{label} (wcEcoli)',
            linestyle='--', linewidth=1, color=color)
    ax.legend()


def get_submass_comparison_plot(
        data, colors, vivarium_agent='', wcecoli_agent=''):
    fig, ax = plt.subplots(figsize=(8, 11))
    plot_submass_comparison(
        ax, data, colors, vivarium_agent, wcecoli_agent)
    return fig
