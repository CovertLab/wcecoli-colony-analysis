from typing import Dict, Sequence, Union, List, Tuple, cast

from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import Colormap
import numpy as np


SerializedField = Union[
    Sequence[Sequence[int]], Sequence[Sequence[float]]]
DEFAULT_COLORMAP = cm.get_cmap('Greys')  # type: ignore
# This avoids plotting data in white
MIN_COLOR_NORMALIZED = 0.2


def get_enviro_sections_plot(
        fields_ts_list: List[Dict[int, Dict[str, SerializedField]]],
        bounds: Sequence,
        section_location: float = 0.5,
        cmap: Colormap = DEFAULT_COLORMAP,
        fontsize: float = 36,
        ) -> Tuple[plt.Figure, dict]:
    '''Get Environment Cross-Sections Plot

    Args:
        fields_ts_list: List of fields timeseries, one timeseries per
            replicate.
        bounds: The (x, y) bounds of the environment.
        section_location: Location along vertical axis of cross-section
            slice. Expressed as a fraction of the vertical height of the
            environment.
        cmap: Colormap from which to draw the colors used to represent
            each timepoint.
        fontsize: Size of all text in figure.

    Returns:
        Tuple of the generated figure and a dictionary of statistics.
        The dictionary maps from times to tuples of the first and third
        quartiles of data across replicates.
    '''
    sorted_times = sorted(fields_ts_list[0].keys())
    some_timepoint = fields_ts_list[0][sorted_times[0]]
    some_field = some_timepoint[list(some_timepoint.keys())[0]]
    num_bins = np.array(some_field).shape  # type: ignore
    bin_width = bounds[0] / num_bins[0]
    x = []
    for i in range(num_bins[0]):
        mid = i * bin_width + bin_width / 2
        x.append(mid)

    y_values_ts: Dict[int, Dict[str, List[List[float]]]] = dict()
    all_fields = list(some_timepoint.keys())
    for time in sorted_times:
        y_values = y_values_ts.setdefault(time, dict())
        for field in all_fields:
            y_values[field] = []
            for fields_ts in fields_ts_list:
                matrix = np.array(  # type: ignore
                    fields_ts[time][field])
                # Skip infinite fields
                if np.any(np.isinf(matrix)):  # type: ignore
                    continue
                assert matrix.shape == num_bins
                section_index = int(section_location * matrix.shape[0])
                section = matrix[section_index, :]
                y_values[field].append(section.tolist())
    field_names = list(y_values_ts[sorted_times[0]].keys())
    num_fields = len(field_names)
    figsize=(num_bins[1] * 0.8, num_fields * 4)
    fig, tmp_axes = plt.subplots(
        num_fields, sharex=True, figsize=figsize)
    # We assume axes is a list
    if num_fields == 1:
        axes: List[plt.Axes] = [cast(plt.Axes, tmp_axes)]
    else:
        axes = tmp_axes
    stats = {}
    for time in sorted_times:
        color_normalized = (
            (1 - MIN_COLOR_NORMALIZED)
            * (time / sorted_times[-1])
            + MIN_COLOR_NORMALIZED)
        color = cmap(color_normalized)
        y_values = y_values_ts[time]
        for field_i, (field, y_list) in enumerate(y_values.items()):
            y_matrix = np.array(y_list)
            median = np.median(y_matrix, axis=0)
            q25, q75 = np.percentile(y_matrix, [25, 75], axis=0)
            stats[time] = q25, median, q75
            ax = axes[field_i]
            ax.plot(  # type: ignore
                x, median, 'o', color=color, linestyle='-',
                label='{:.0f}s'.format(time))
            ax.fill_between(  # type: ignore
                x, q25, q75, color=color, alpha=0.2, edgecolor='none')
    for i, ax in enumerate(axes):
        ax.tick_params(  # type: ignore
            axis='both', which='major', labelsize=fontsize)
        ax.set_xticks(np.linspace(0, bounds[0], 6))
        if num_fields != 1:
            ax.set_title(  # type: ignore
                field_names[i], fontsize=fontsize)
        if i == num_fields - 1:
            ax.set_xlabel(  # type: ignore
                'Environment Horizontal Axis ($\\mu m$)',
                fontsize=fontsize)
            ax.legend(  # type: ignore
                bbox_to_anchor=(1.05, 0.5), loc='center left',
                prop={'size': fontsize})
    # Make y label centered across all subplots
    super_ax = fig.add_subplot(  # type: ignore
        111, xticks=[], yticks=[], frameon=False)
    super_ax.set_ylabel(  # type: ignore
        'Concentration ($mM$)', labelpad=75, fontsize=fontsize)
    fig.tight_layout()
    return fig, stats
