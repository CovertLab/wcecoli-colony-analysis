import argparse
import json
from typing import Callable, List

import numpy as np
from scipy import stats as scipy_stats
from scipy.constants import N_A


# From the "expression_distributions" stats in figs42. In counts/fL
EXPRESSION_IQRS = {
    'AmpC': (68.36393149003717, 85.45342670347938, 99.41098877494602),
    'AcrAB-TolC': (
        18.06286132322007, 16.957394256525866, 19.251988781503236),
}


def analyze_expression_distributions_stats(stats: dict) -> dict:
    summary = {}
    for protein, (minimum, q1, q2, q3, maximum) in stats.items():
        summary[protein] = {
            'Median Concentration (counts/fL)': q2,
            'IQR': (np.array(q3) - np.array(q1)).tolist(),
            'Minimum': minimum,
            'Maximum': maximum,
        }
    return summary


def analyze_growth_fig_stats(stats: dict) -> dict:
    summary = {}
    for condition, (q1, q2, q3) in stats.items():
        summary[condition] = {
            'Median Total Mass Fold-Change': q2[-1] / q2[0],
            'Maximum IQR / Median': (
                (np.array(q3) - np.array(q1)) / q2).max(),
        }
    return summary


def analyze_threshold_scan_stats(stats: dict) -> dict:
    summary = analyze_growth_fig_stats(stats)
    for threshold, (q1, q2, q3) in stats.items():
        end_summary = {
            'Median Final Mass (fg)': q2[-1],
            'Final Mass IQR (fg)': q3[-1] - q1[-1],
        }
        summary[threshold].update(end_summary)
    return summary


def analyze_enviro_section_stats(stats: dict) -> dict:
    summary = {}
    for timepoint, (q1, q2, q3) in stats.items():
        summary[timepoint] = {
            'Maximum IQR (mM)': (np.array(q3) - np.array(q1)).max(),
            'Maximum Median (mM)': max(q2),
            'Minimum Median (mM)': min(q2),
        }
    return summary


def _u_power_centrality(
        num_a: int, num_b: int, colony_radius: float,
        get_prob_a: Callable[[float], float], iters: int = 10000,
        alpha: float = 0.05, seed: int = 530) -> float:
    num_points = num_a + num_b
    p_values = []
    random = np.random.default_rng(seed)  # type: ignore
    for _ in range(iters):
        points = random.uniform(
            -colony_radius, colony_radius, size=(10 * num_points, 2))
        mask = np.linalg.norm(  # type: ignore
            points, ord=2, axis=1) <= colony_radius
        points = points[mask]

        a_points: List[np.ndarray] = []
        b_points: List[np.ndarray] = []
        for point in points:
            dist = np.linalg.norm(point, ord=2)  # type: ignore
            prob_a = get_prob_a(dist)
            is_a = random.random() < prob_a  # type: ignore
            if is_a and len(a_points) < num_a:
                a_points.append(point)
            elif len(b_points) < num_b:
                b_points.append(point)

        a_dists = np.linalg.norm(  # type: ignore
            a_points, ord=2, axis=1)
        b_dists = np.linalg.norm(  # type: ignore
            b_points, ord=2, axis=1)

        _, p_value = scipy_stats.mannwhitneyu(
            a_dists,
            b_dists,
            alternative='two-sided',
        )
        p_values.append(p_value)

    p_arr = np.array(p_values)
    p = (p_arr < alpha).sum() / len(p_arr)
    return p


def analyze_centrality_stats(stats: dict) -> dict:
    summary = {}
    survive_q1, survive_q2, survive_q3 = np.percentile(
        stats['survive_distances'], [25, 50, 75])
    die_q1, die_q2, die_q3 = np.percentile(
        stats['die_distances'], [25, 50, 75])
    summary['survive'] = {
        'Median Euclidian distance from center (um)': survive_q2,
        'IQR (um)': survive_q3 - survive_q1,
    }
    summary['die'] = {
        'Median Euclidian distance from center (um)': die_q2,
        'IQR (um)': die_q3 - die_q1,
    }
    u_stat, p_value = scipy_stats.mannwhitneyu(
        stats['survive_distances'],
        stats['die_distances'],
        alternative='two-sided',
    )
    summary['hypothesis testing'] = {
        'Mann-Whitney U statistic': u_stat,
        'Mann-Whitney p-value': p_value,
        'Power (alpha=0.2)for 0.5 diff in p(death)': _u_power_centrality(
            len(stats['survive_distances']),
            len(stats['die_distances']),
            10,
            lambda x: 0.25 + x / 10 / 2,
            alpha=0.2,
        ),
    }
    return summary


def analyze_growth_snapshot_stats(stats: dict) -> dict:
    summary = {}
    for condition, condition_stats in stats.items():
        final_agent_counts = []
        for replicate_stats in condition_stats.values():
            times = [
                float(time)
                for time in replicate_stats['agents'].keys()]
            initial_agents = replicate_stats['agents'][str(min(times))]
            final_agents = replicate_stats['agents'][str(max(times))]
            assert initial_agents == 1
            final_agent_counts.append(final_agents)
        assert len(set(final_agent_counts)) == 1
        summary[condition] = {
            'Final number of agents': final_agent_counts[0],
        }
    return summary


def analyze_enviro_heterogeneity_stats(stats: dict) -> dict:
    summary: dict = {}
    replicate_summaries: dict = {}
    for replicate_stats in stats.values():
        for field, field_stats in replicate_stats['fields'].items():
            field_summary = replicate_summaries.setdefault(field, {})
            times = [float(time) for time in field_stats]
            for name, func in {'initial': min, 'final': max}.items():
                field_summary.setdefault(
                    name,
                    {
                        key: []
                        for key in ('min', 'median', 'max', 'iqr')
                    },
                )
                time = func(times)
                time_min, q1, q2, q3, time_max = field_stats[str(time)]
                field_summary[name]['min'].append(time_min)
                field_summary[name]['median'].append(q2)
                field_summary[name]['max'].append(time_max)
                field_summary[name]['iqr'].append(q3 - q1)

    for field, field_summary in replicate_summaries.items():
        summary[field] = {}
        for time, time_summary in field_summary.items():
            summary[field][time] = {}
            for key, array in time_summary.items():
                q1, q2, q3 = np.percentile(array, [25, 50, 75])
                summary[field][time][key] = {
                    'Median across replicates (mM)': q2,
                    'IQR across replicates (mM)': q3 - q1,
                }
    return summary


def _u_power_concentrations(
        num_a: int, num_b: int, a_stdev: float, b_stdev: float,
        diff: float, iters: int = 10000, seed: int = 620,
        alpha: float = 0.05) -> float:
    random = np.random.default_rng(seed)  # type: ignore
    a_values = random.normal(size=(iters, num_a), scale=a_stdev)
    b_values = random.normal(
        size=(iters, num_b), loc=diff, scale=b_stdev)
    p_values = []
    for i in range(iters):
        a = a_values[i, :]
        b = b_values[i, :]
        _, p_value = scipy_stats.mannwhitneyu(a, b, alternative='two-sided')
        p_values.append(p_value)
    p_arr = np.array(p_values)
    power = (p_arr < alpha).sum() / len(p_arr)
    return power


def analyze_dotplot_stats(stats: dict) -> dict:
    summary: dict = {}
    stdevs = {
        # Convert IQRs from counts/fL to mM
        protein: (np.array(iqrs) * 1e15 / 1e3 / N_A).mean() / (  # type: ignore
            scipy_stats.norm.ppf(0.75) - scipy_stats.norm.ppf(0.25))
        for protein, iqrs in EXPRESSION_IQRS.items()
    }
    for protein, protein_stats in stats.items():
        summary[protein] = {}
        for status, concentrations in protein_stats.items():
            q1, q2, q3 = np.percentile(concentrations, [25, 50, 75])
            summary[protein][status] = {
                'Median (mM)': q2,
                'IQR (mM)': q3 - q1,
            }
        u_stat, p_value = scipy_stats.mannwhitneyu(
            protein_stats['live'],
            protein_stats['dead'],
            alternative='two-sided',
        )
        summary[protein]['Mann-Whitney U-test'] = {
            'U-statistic': u_stat,
            'p-value': p_value,
            'power': _u_power_concentrations(
                len(protein_stats['live']),
                len(protein_stats['dead']),
                stdevs[protein], stdevs[protein], 1e-4,
            ),
        }
    return summary


SECTION_ANALYZER_MAP = {
    'expression_distributions': analyze_expression_distributions_stats,
    'growth_fig': analyze_growth_fig_stats,
    'threshold_scan': analyze_threshold_scan_stats,
    'enviro_section': analyze_enviro_section_stats,
    'centrality': analyze_centrality_stats,
    'growth_snapshots': analyze_growth_snapshot_stats,
    'enviro_heterogeneity': analyze_enviro_heterogeneity_stats,
    'dotplots': analyze_dotplot_stats,
}


def main(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'stats_json',
        help='Path to stats JSON file')
    parser.add_argument(
        '-o', '--out',
        default='summary_stats.json',
        help='Path to write summary stats to')

    args = parser.parse_args(args)

    with open(args.stats_json, 'r') as f:
        stats = json.load(f)

    summary = {}
    for section, analyzer in SECTION_ANALYZER_MAP.items():
        summary[section] = analyzer(stats[section])

    with open(args.out, 'w') as f:
        json.dump(summary, f, indent=4)


if __name__ == '__main__':
    main()
