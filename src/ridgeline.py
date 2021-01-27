from typing import (
    Dict, Sequence, Union, Tuple, List, Optional, Any, Iterable, cast)

from scipy.stats.kde import gaussian_kde
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes


Color = Union[str, Tuple[float, float, float, float]]
Number = Union[float, int, np.float64]


def flatten(value: Any) -> list:
    '''Flatten any iterable into a single list of its base elements.

    Strings are considered base elements, not iterables, for this
    function. If ``value`` is not iterable or a string, a list is
    returned containing only ``value``. Ordering within each iterable is
    maintained in the final list, and flattening occurs in a depth-first
    fashion. For example:

    >>> flatten([[1, 2], [3, [4, 5]]])
    [1, 2, 3, 4, 5]

    Args:
        value: The value to flatten.
    Returns:
        A flat list containing all the base elements of ``value``.
    '''
    try:
        iter(value)
    except TypeError:
        return [value]
    if isinstance(value, str):
        return [value]
    flat = []
    lst = cast(Iterable, value)
    for elem in lst:
        flat.extend(flatten(elem))
    return flat


def get_ridgeline_plot(
        replicates: Iterable[Tuple[
            Dict[str, Sequence[Number]],
            Color,
        ]],
        point_alpha: float = 1,
        num_bins: int = 20,
        overlap: float = 0.2,
        horizontal_extra: float = 0.2,
        jitter: Optional[float] = None,
        figsize: Optional[Tuple[Number, Number]] = None,
        x_label: str = '',
        y_label: str = '',
        ) -> plt.Figure:
    '''Generate a ridgeline plot.

    This convenience function generates a figure and uses
    :py:func:`plot_ridgeline` to create a ridgeline plot on the figure.
    Only the arguments specific to this function are documented here.
    For the rest, see :py:func:`plot_ridgeline`.

    Args:
        replicates: An iterable of tuples. Each tuple contains first a
            data dictionary as described in :py:func:`plot_ridgeline`
            and second a color. Each tuple describes the data from a
            replicate. All replicates will be plotted together to show
            variation.
        figsize: Tuple indicating the figure size in inches as ``(width,
        height)``.
        x_label: Label for x axis.
        y_label: Label for y axis.

    Returns:
        A figure with the ridgeline plot.
    '''
    if figsize is None:
        fig, ax = plt.subplots()
    else:
        fig, ax = plt.subplots(figsize=figsize)
    data = [replicate[0] for replicate in replicates]
    colors = [replicate[1] for replicate in replicates]
    plot_ridgeline(
        data, ax, colors, 0.2, point_alpha, num_bins, overlap,
        horizontal_extra, jitter)
    if x_label:
        ax.set_xlabel(x_label)
    if y_label:
        ax.set_ylabel(y_label)
    for spine_name in ('top', 'right'):
        ax.spines[spine_name].set_visible(False)  # type: ignore
    fig.tight_layout()
    return fig


def _calculate_density_curves(
        data: Sequence[Dict[str, Sequence[Number]]],
        x_values: Union[Sequence[Number], np.ndarray],
        overlap: float):
    y_values: Dict[str, List[float]] = {}
    density_curves: Dict[str, List[Sequence[float]]] = {}
    # Reverse order since we plot from bottom to top
    for data_dict in data:
        y = 0.
        for y_label in reversed(data_dict.keys()):
            data_values = data_dict[y_label]
            if len(data_values) > 1:
                pdf = gaussian_kde(data_values)
                density_curve = pdf(x_values)
            elif len(data_values) == 1:
                density_curve = np.zeros(len(x_values))
                for i, x_value in enumerate(x_values):
                    if x_value >= data_values[0]:
                        density_curve[i] = 1
                        break
            else:
                density_curve = np.zeros(len(x_values))
            density_curves.setdefault(y_label, []).append(density_curve)
            y_values.setdefault(y_label, []).append(y)
            y += max(density_curve) * (1 - overlap)
    max_y_values = {
        y_label: max(values)
        for y_label, values in y_values.items()
    }
    return max_y_values, density_curves



def plot_ridgeline(
        data: Sequence[Dict[str, Sequence[Number]]],
        ax: Axes,
        colors: Sequence[Color],
        fill_alpha: float = 0.2,
        point_alpha: float = 1,
        num_bins: int = 20,
        overlap: float = 0.2,
        horizontal_extra: float = 0.2,
        jitter: Optional[float] = None,
        data_bounds: Optional[Tuple[Number, Number]] = None,
        ) -> None:
    '''Plot data as a ridgeline plot.

    Args:
        data: List of mappings from strings to sequences of values. One
            distribution will be plotted for each sequence of values and
            labeled with the corresponding string. Each mapping
            represents a replicate.
        ax: Axes on which to plot.
        colors: The colors to plot each replicate in.
        fill_alpha: Alpha value to control transparency of fill color.
        point_alpha: Alpha value to control transparency of point color.
        num_bins: Number of equally-sized bins into which the values
            will be grouped. The averages of the values in each bin
            determine the curve height at the middle of the bin on the
            x-axis. This is equivalent to the bins on a histogram.
        overlap: The fraction by which distributions are allowed to
            overlap. For example, with an overlap of ``0.2``, the first
            distribution would be at :math:`y=0` and the next
            distribution would be at :math:`y=0.8`. Note that since not
            all distributions will reach a height of :math:`1`, setting
            a nonzero overlap does not guarantee that any distributions
            will actually overlap.
        horizontal_extra: The fraction of the range of ``data`` which
            should be added to the negative and positive x axis. This
            adds extra space to the left and right of distributions.
        jitter: To let the viewer distinguish between multiple nearby
            points, we randomly add small values to point positions.
            ``jitter`` is the maximum absolute value of these
            perturbations, which are uniformly chosen from the range
            :math:`(-jitter, jitter)`. If ``None``, ``jitter`` will be
            set to one tenth of the maximum value in ``data``. To
            disable jitter, set ``jitter=0``.
        data_bounds: The minimum and maximum values present in the data.
            If ``None``, the range is calculated from ``data``. This may
            be useful when ``data`` is only some of the data that will
            be plotted, e.g. with replicates.
    '''
    if data_bounds is None:
        flat_data = flatten([data_elem.values() for data_elem in data])
        data_min = min(flat_data)
        data_max = max(flat_data)
    else:
        data_min, data_max = data_bounds
    data_range = data_max - data_min
    if jitter is None:
        jitter = data_range / 10
    extra = data_range * horizontal_extra
    x_values = np.linspace(
        data_min - extra, data_max + extra, num_bins)
    y_values, density_curves = _calculate_density_curves(
        data, x_values, overlap)
    num_replicates = len(density_curves[list(y_values.keys())[0]])
    assert len(colors) == num_replicates

    for i, y_label in enumerate(y_values):
        y = y_values[y_label]
        assert len(density_curves[y_label]) == num_replicates
        for j, density_curve in enumerate(density_curves[y_label]):
            color = colors[j]
            zorder = num_replicates - i + 1
            # Below, we cast colors to strings to satisfy mypy, which
            # doesn't understand using tuples as colors.
            ax.plot(
                x_values, density_curve + y,
                color=cast(str, color), linestyle='-', zorder=zorder)
            # Mypy incorrectly complains that Axes has no fill_between
            # attribute.
            ax.fill_between(  # type: ignore
                x_values, np.ones(num_bins) * y, density_curve + y,
                color=cast(str, color), zorder=-1, alpha=fill_alpha)
            points = np.array(
                data[j][y_label], dtype=float)  # type: ignore
            # We disable pylint's no-member check because pylint doesn't
            # recognize that np.random.random is valid.
            # pylint: disable=no-member
            points += (
                np.random.random(len(points)) - 0.5  # type: ignore
            ) * 2 * jitter
            # pylint: enable=no-member
            ax.scatter(points,  # type: ignore
                np.ones(len(points)) * y, color=cast(str, color),
                marker='|', s=20, linewidths=0.1, zorder=zorder,
                alpha=point_alpha)
    ax.set_yticks(list(y_values.values()))
    ax.set_yticklabels(list(y_values.keys()))
