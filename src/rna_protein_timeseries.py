from matplotlib import pyplot as plt
from vivarium.library.topology import get_in


def plot_rna_protein_timeseries(
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
        ):
    # data: The raw simulation data for a single agent.
    # name_path_map: Map from the name to show on the plot to the path
    # in data.
    # names_layout: Nested lists forming a rectangle where each cell is
    # either None (for no plot) or a name from name_path_map. Each list
    # represents a row.

    n_cols = len(names_layout[0])
    n_rows = len(names_layout)
    fig = plt.figure(figsize=(n_cols * col_width, n_rows * row_height))
    grid = plt.GridSpec(
        ncols=n_cols, nrows=n_rows, wspace=0.4, hspace=1.5)
    time = list(data.keys())
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
            ax.set_xlim([time[0], time[-1]])
            ax.xaxis.get_offset_text().set_fontsize(tick_label_size)
            ax.yaxis.get_offset_text().set_fontsize(tick_label_size)
            if i == n_rows - 1 or names_layout[i + 1][j] is None:
                ax.set_xlabel(x_label, fontsize=axis_label_size)

            paths = name_paths_map[name]
            series = []
            for timepoint in data.values():
                value = 0
                for (path, multiplier) in paths:
                    agent_data = get_in(timepoint, agent_path)
                    value += get_in(agent_data, path) * multiplier
                series.append(value)
            print(name, series)
            ax.plot(time, series, linewidth=linewidth)
    fig.tight_layout()
    return fig