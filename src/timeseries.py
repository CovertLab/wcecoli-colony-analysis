from matplotlib import pyplot as plt
import numpy as np
from vivarium.library.topology import get_in
from vivarium.plots.simulation_output import set_axes


def get_timeseries_plot(
        data,
        name_paths_map,
        name_y_label_map,
        names_layout,
        col_width=4,
        row_height=2,
        tick_label_size=12,
        title_size=16,
        axis_label_size=12,
        x_label='time (s)',
        agent_path=tuple(),
        linewidth=3,
        horizontal_lines=None,
        min_value=-np.inf,
        ):
    '''Create one or more timeseries plots.

    Args:
        data: Raw data dictionary.
        name_paths_map: Mapping from a name of a variable to the path at
            which that variable can be found in data. If ``agent_path``
            is specified, then paths are prepended with ``agent_path``
            first. The name will be shown as the title of the plot.
        name_y_label_map: Map from variable names to the y-axis label to
            show for the variable. If a variable is not present, it no
            y-axis label will be shown for its plot.
        names_layout: Nested lists forming a matrix where each cell is
            either None (for no plot) or a name from ``name_path_map``.
            Each list represents a row in the final layout.
        col_width: Width of a column of plots.
        row_height: Height of a row of plots.
        tick_label_size: Size of tick labels.
        title_size: Font size for titles.
        axis_label_size: Font size for x- and y-axis labels.
        x_label: Label for the x-axis.
        agent_path: Path at which the agent to plot data for can be
            found in ``data``.
        horizontal_lines: Optional dictionary specifying horizontal
            lines to add to the plot. The dictionary maps from variable
            names to a tuple of line specifiers. Each line specifier is
            a 3-tuple of (label, y-axis-location, color). The label will
            be show in the legend to identify the line.
        min_value: If specified, all values smaller than ``min_value``
            will be excluded from the plot. All plots must have the same
            times after exclusion.

    Returns:
        The plot.
    '''
    horizontal_lines = horizontal_lines or {}
    n_cols = len(names_layout[0])
    n_rows = len(names_layout)
    fig = plt.figure(figsize=(n_cols * col_width, n_rows * row_height))
    grid = plt.GridSpec(
        ncols=n_cols, nrows=n_rows, wspace=0.4, hspace=1.5)
    times = []
    for i in range(n_rows):
        for j in range(n_cols):
            name = names_layout[i][j]
            if name is None:
                continue
            ax = fig.add_subplot(grid[i, j])
            ax.set_title(name)
            ax.title.set_fontsize(title_size)
            for tick_type in ('major', 'minor'):
                ax.tick_params(
                    axis='both',
                    which=tick_type,
                    labelsize=tick_label_size,
                )
            y_label = name_y_label_map[name]
            if y_label:
                ax.set_ylabel(y_label, fontsize=axis_label_size)
            ax.xaxis.get_offset_text().set_fontsize(tick_label_size)
            ax.yaxis.get_offset_text().set_fontsize(tick_label_size)
            if i == n_rows - 1 or names_layout[i + 1][j] is None:
                ax.set_xlabel(x_label, fontsize=axis_label_size)

            paths = name_paths_map[name]
            series = []
            series_times = []
            for time, timepoint in data.items():
                agent_data = get_in(timepoint, agent_path)
                if not agent_data:
                    continue
                value = 0
                for (path, multiplier) in paths:
                    value += get_in(agent_data, path) * multiplier
                if value > min_value:
                    series.append(value)
                    series_times.append(time)
            if not times:
                times = series_times
            else:
                assert times == series_times
            ax.set_xlim([times[0], times[-1]])
            ax.plot(times, series, linewidth=linewidth)
            for label, y_loc, color in horizontal_lines.get(
                    name, tuple()):
                ax.hlines(
                    y_loc, times[0], times[-1], label=label,
                    linestyles='dashed', colors=[color])
                ax.legend()
            set_axes(ax, True, sci_notation=True)
    fig.tight_layout()
    return fig
