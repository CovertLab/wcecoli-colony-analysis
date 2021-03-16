import argparse
import json
from typing import Callable, List

import numpy as np
from scipy import stats as scipy_stats


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
            'Maximum IQR (fg)': (np.array(q3) - np.array(q1)).max(),
        }
    return summary


def analyze_enviro_section_stats(stats: dict) -> dict:
    summary = {}
    for timepoint, (q1, q2, q3) in stats.items():
        summary[timepoint] = {
            'Maximum IQR (mM)': (np.array(q3) - np.array(q1)).max(),
            'Maximum Median (mM)': q2.max(),
            'Minimum Median (mM)': q2.min(),
        }
    return summary


def analyze_threshold_scan_stats(stats: dict) -> dict:
    summary = {}
    for threshold, (q1, q2, q3) in stats.items():
        summary[threshold] = {
            'Median Final Mass (fg)': q2[-1],
            'Final Mass IQR (fg)': (np.array(q3) - np.array(q1))[-1],
        }
    return summary


def _mann_whitney_u_power(
        num_a: int, num_b: int, colony_radius: float,
        get_prob_a: Callable[[float], float], iters: int = 10000,
        alpha: float = 0.05) -> float:
    num_points = num_a + num_b
    p_values = []
    for _ in range(iters):
        points = np.random.uniform(
            -colony_radius, colony_radius, size=(2, 10 * num_points))
        mask = np.linalg.norm(  # type: ignore
            points, ord=2, axis=1) <= colony_radius
        points = points[mask]

        a_points: List[np.ndarray] = []
        b_points: List[np.ndarray] = []
        for point in points:
            dist = np.linalg.norm(point, ord=2)  # type: ignore
            prob_a = get_prob_a(dist)
            is_a = np.random.random() < prob_a  # type: ignore
            if is_a and len(a_points) < num_a:
                a_points.append(point)
            elif len(b_points) < num_b:
                b_points.append(point)

        a_dists = np.linalg.norm(  # type: ignore
            a_points, ord=2, axis=1)
        b_dists = np.linalg.norm(  # type: ignore
            b_points, ord=2, axis=1)

        _, p_value = scipy_stats.mannwhitneyu(a_dists, b_dists)
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
        'Median Euclidian distance from center': survive_q2,
        'IQR': survive_q3 - survive_q1,
    }
    summary['die'] = {
        'Median Euclidian distance from center': die_q2,
        'IQR': die_q3 - die_q1,
    }
    u_stat, p_value = scipy_stats.mannwhitneyu(
        stats['survive_distances'], stats['die_distances'])
    summary['hypothesis testing'] = {
        'Mann-Whitney U statistic': u_stat,
        'Mann-Whitney p-value': p_value,
        'Power for 0.01 diff in p(death)': _mann_whitney_u_power(
            len(stats['survive_distances']),
            len(stats['die_distances']),
            10,
            lambda x: 0.495 + x / 10 / 100
        ),
    }
    return summary


SECTION_ANALYZER_MAP = {
    'expression_distributions': analyze_expression_distributions_stats,
    'growth_fig': analyze_growth_fig_stats,
    'enviro_section': analyze_enviro_section_stats,
    'centrality': analyze_centrality_stats,
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
