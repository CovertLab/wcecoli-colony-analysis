import argparse
import json

import numpy as np


def analyze_expression_distributions_stats(stats: dict) -> dict:
    summary = {}
    for protein, (minimum, q1, q2, q3, maximum) in stats.items():
        summary[protein] = {
            'Median Concentration (counts/fL)': q2,
            'IQR': q3 - q1,
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
    for timepoint, (q1, q3) in stats.items():
        summary[timepoint] = {
            'Maximum IQR (mM)': (np.array(q3) - np.array(q1)).max(),
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


SECTION_ANALYZER_MAP = {
    'expression_distributions': analyze_expression_distributions_stats,
    'growth_fig': analyze_growth_fig_stats,
    'enviro_section': analyze_enviro_section_stats,
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
