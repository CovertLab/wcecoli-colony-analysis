from typing import Dict, Sequence, Union, List

from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import Colormap
import numpy as np


SerializedField = Union[
    Sequence[Sequence[int]], Sequence[Sequence[float]]]
DEFAULT_COLORMAP = cm.get_cmap('Greys')
# This avoids plotting data in white
MIN_COLOR_NORMALIZED = 0.2


def get_enviro_sections_plot(
        fields_ts: Dict[int, Dict[str, SerializedField]],
        bounds: Sequence,
        section_location: float = 0.5,
        flat_bins: bool = False,
        cmap: Colormap = DEFAULT_COLORMAP) -> plt.Figure:
    sorted_times = sorted(fields_ts.keys())
    some_timepoint = fields_ts[sorted_times[0]]
    some_field = some_timepoint[list(some_timepoint.keys())[0]]
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

    y_values_ts: Dict[int, Dict[str, List[float]]] = dict()
    for time in sorted_times:
        fields = fields_ts[time]
        y_values = y_values_ts.setdefault(time, dict())
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
    field_names = list(y_values_ts[sorted_times[0]].keys())
    num_fields = len(field_names)
    figsize=(num_bins[1], num_fields * 4)
    fig, axes = plt.subplots(
        num_fields, sharex=True, figsize=figsize)
    if num_fields == 1:
        # We assume axes is a list
        axes = [axes]
    for time in sorted_times:
        color_normalized = (
            (1 - MIN_COLOR_NORMALIZED)
            * (time / sorted_times[-1])
            + MIN_COLOR_NORMALIZED)
        color = cmap(color_normalized)
        y_values = y_values_ts[time]
        for field_i, (field, y) in enumerate(y_values.items()):
            ax = axes[field_i]
            if flat_bins:
                ax.plot(x, y, color=color, label='{}s'.format(time))
            else:
                ax.plot(x, y, 'o', color=color, linestyle='-',
                        label='{}s'.format(time))
    for i, ax in enumerate(axes):
        if num_fields != 1:
            ax.set_title(field_names[i])
        if i == num_fields - 1:
            ax.set_xlabel('Environment Horizontal Axis ($\\mu m$)')
            ax.legend(bbox_to_anchor=(1.5, 1), loc='upper right')
    # Make y label centered across all subplots
    super_ax = fig.add_subplot(111, xticks=[], yticks=[], frameon=False)
    super_ax.set_ylabel('Concentration ($mM$)', labelpad=75)
    fig.tight_layout()
    return fig
