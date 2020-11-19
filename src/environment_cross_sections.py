from typing import Dict, Sequence, Union

from matplotlib import pyplot as plt
import numpy as np


SerializedField = Union[
    Sequence[Sequence[int]], Sequence[Sequence[float]]]


def get_enviro_sections_plot(
        fields: Dict[str, SerializedField],
        bounds: Sequence,
        section_location: float = 0.5,
        flat_bins: bool = False) -> plt.Figure:
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
    y_values = {}
    for field, field_vals in fields.items():
        matrix = np.array(field_vals)
        # Skip infinite fields
        if np.any(np.isinf(matrix)):
            continue
        assert matrix.shape == num_bins
        section_index = int(section_location * matrix.shape[0])
        section = matrix[section_index, :]
        # Value for start and end points
        y = []
        if flat_bins:
            for val in section:
                y += [val, val]
        else:
            y = section.tolist()
        y_values[field] = y
    figsize=(num_bins[1], len(y_values))
    fig, axes = plt.subplots(
        nrows=len(y_values), sharex=True, figsize=figsize)
    for i, (field, y) in enumerate(y_values.items()):
        ax = axes[i]
        if flat_bins:
            ax.plot(x, y)
        else:
            ax.plot(x, y, 'o', color='black', linestyle='-')
        ax.set_title(field)
        if i == len(y_values) - 1:
            ax.set_xlabel('Environment Horizontal Axis ($\\mu m$)')
    # Make y label centered across all subplots
    super_ax = fig.add_subplot(111, xticks=[], yticks=[], frameon=False)
    super_ax.set_ylabel('Concentration ($mM$)', labelpad=75)
    fig.tight_layout()
    return fig
