import itertools
from typing import Dict, Sequence, Union
import warnings

from matplotlib import pyplot as plt
import numpy as np


SerializedField = Union[
    Sequence[Sequence[int]], Sequence[Sequence[float]]]


LINE_STYLES = ('-', ':', '-.', '--')
# Generated using https://medialab.github.io/iwanthue/
COLORS = (
    '#cd50be', '#70d658', '#7544cc', '#cdd655', '#512770', '#d09533',
    '#727bce', '#d74830', '#64d0aa', '#ce4974', '#63893b', '#936088',
    '#bdd3ad', '#39344c', '#c59f78', '#72abc2', '#b16241', '#d3a6cc',
    '#425439', '#682f2d',
)
STYLES = tuple(itertools.product(LINE_STYLES, COLORS))


def _get_shape(matrix: SerializedField):
    num_rows = len(matrix)
    num_cols = len(matrix[0])
    return num_rows, num_cols


def get_enviro_sections_plot(
    fields: Dict[str, SerializedField],
    bounds: Sequence,
    section_location: float = 0.5,
    flat_bins: bool = False,
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(30, 30))
    plot_enviro_sections(ax, fields, bounds, section_location, flat_bins)
    fig.tight_layout()
    return fig


def plot_enviro_sections(
    ax: plt.Axes,
    fields: Dict[str, SerializedField],
    bounds: Sequence,
    section_location: float = 0.5,
    flat_bins: bool = False,
) -> None:
    some_field = fields[list(fields.keys())[0]]
    num_bins = np.array(some_field).shape
    bin_width = bounds[0] / num_bins[0]
    x = []
    if flat_bins:
        for i in range(num_bins[0]):
            start = i * bin_width
            end = (i + 1) * bin_width
            x += [start, end]
    else:
        for i in range(num_bins[0]):
            mid = i * bin_width + bin_width / 2
            x.append(mid)
    i = 0
    for field, field_vals in fields.items():
        matrix = np.array(field_vals)
        # Skip infinite fields
        if np.any(np.isinf(matrix)):
            continue
        assert matrix.shape == num_bins
        section_index = int(section_location * matrix.shape[0])
        section = matrix[section_index, :]
        normalized = section - section[0]
        # Value for start and end points
        y = []
        if flat_bins:
            for val in normalized:
                y += [val, val]
        else:
            y = normalized.tolist()
        linestyle, color = STYLES[i % len(STYLES)]
        if i > len(STYLES) - 1:
            warnings.warn(
                'More than {} non-inf fields so reusing style'.format(
                    len(STYLES))
            )
        if flat_bins:
            ax.plot(x, y, label=field, linestyle=linestyle, color=color)
        else:
            ax.plot(x, y, 'o', label=field, linestyle=linestyle,
                    color=color)
        i += 1
    ax.legend()
    ax.set_xlabel('Environment Horizontal Axis ($\mu m$)')
    ax.set_ylabel('Change in Concentration from $x = 0$ ($mM$)')
