from typing import (
    Dict, Sequence, Union, Tuple, List, Optional, Any, Iterable, cast)

from scipy.stats.kde import gaussian_kde
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes


Color = Union[str, Tuple[float, float, float, float]]
Number = Union[float, int]


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
        num_bins: int = 20,
        overlap: float = 0.2,
        horizontal_extra: float = 0.2,
        jitter: Optional[float] = None,
        figsize: Optional[Tuple[Number, Number]] = None,
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

    Returns:
        A figure with the ridgeline plot.
    '''
    if figsize is None:
        fig, ax = plt.subplots()
    else:
        fig, ax = plt.subplots(figsize=figsize)
    flat_data: List[Number] = []
    for data, _ in replicates:
        flat_data.extend(flatten(data.values()))
    data_bounds = min(flat_data), max(flat_data)
    for data, color in replicates:
        plot_ridgeline(
            data, ax, color, 0.2, color, color, num_bins, overlap,
            horizontal_extra, jitter, data_bounds)
    return fig


def plot_ridgeline(
        data: Dict[str, Sequence[Union[Number]]],
        ax: Axes,
        fill_color: Optional[Color] = 'gray',
        fill_alpha: float = 0.2,
        line_color: Optional[Color] = 'black',
        point_color: Optional[Color] = 'black',
        num_bins: int = 20,
        overlap: float = 0.2,
        horizontal_extra: float = 0.2,
        jitter: Optional[float] = None,
        data_bounds: Optional[Tuple[Number, Number]] = None,
        ) -> None:
    '''Plot data as a ridgeline plot.

    Args:
        data: A mapping from strings to sequences of values. One
            distribution will be plotted for each sequence of values and
            labeled with the corresponding string.
        ax: Axes on which to plot.
        fill_color: Color to show under each distribution. If ``None``,
            the area under each distribution will not be filled.
        fill_alpha: Alpha value to control transparency of fill color.
        line_color: Color of curve that forms the top of each
            distribution. If ``None``, no curve is shown.
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
        point_color: The color of markers along the bottom of each
            distribution showing the actual point values. If ``None``,
            no markers are plotted.
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
    # Mypy doesn't recognize that a view of dictionary values is valid
    # input for np.array()
    if data_bounds is None:
        flat_data = flatten(data.values())
        data_min = min(flat_data)
        data_max = max(flat_data)
    else:
        data_min, data_max = data_bounds
    if jitter is None:
        jitter = data_max / 10
    data_range = data_max - data_min
    extra = data_range * horizontal_extra
    x_values = np.linspace(
        data_min - extra, data_max + extra, num_bins)
    y_values: List[float] = []
    y_labels: List[str] = []
    for i, (label, values) in enumerate(data.items()):
        pdf = gaussian_kde(values)
        y = i * (1 - overlap)
        y_values.append(y)
        y_labels.append(label)
        density_curve = pdf(x_values)
        zorder = len(data) - i + 1

        # Below, we cast colors to strings to satisfy mypy, which
        # doesn't understand using tuples as colors.
        if line_color is not None:
            ax.plot(
                x_values, density_curve + y,
                color=cast(str, line_color), linestyle='-',
                zorder=zorder)
        if fill_color is not None:
            # Mypy incorrectly complains that Axes has no fill_between
            # attribute.
            ax.fill_between(  # type: ignore
                x_values, np.ones(num_bins) * y, density_curve + y,
                color=cast(str, fill_color), zorder=-1,
                alpha=fill_alpha)
        if point_color is not None:
            point_values = np.array(
                values, dtype=float)  # type: ignore
            # We disable pylint's no-member check because pylint doesn't
            # recognize that np.random.random is valid.
            # pylint: disable=no-member
            point_values += (
                np.random.random(len(values)) - 0.5  # type: ignore
            ) * 2 * jitter
            # pylint: enable=no-member
            ax.scatter(point_values,  # type: ignore
                np.ones(len(values)) * y, color=cast(str, point_color),
                marker='|', zorder=zorder)
    ax.set_yticks(y_values)
    ax.set_yticklabels(y_labels)
