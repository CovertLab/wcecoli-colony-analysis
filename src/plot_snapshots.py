# type: ignore  # pylint: skip-file
'''
This code was adapted from vivarium-cell version 0.0.24. That project is
hosted on GitHub at
https://github.com/vivarium-collective/vivarium-cell/ and is under an MIT
license.

Copyright (c) 2018-2021 The Vivarium Authors, 2021 Christopher Skalnik
'''

from __future__ import absolute_import, division, print_function

import os
import math
import random
import itertools

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
from matplotlib.colors import hsv_to_rgb
from mpl_toolkits.axes_grid1 import (
    make_axes_locatable, anchored_artists)
import numpy as np
from vivarium.library.dict_utils import get_value_from_path
from vivarium.library.topology import get_in


DEFAULT_BOUNDS = [10, 10]

# constants
PI = math.pi

# colors for phylogeny initial agents
HUES = [hue/360 for hue in np.linspace(0,360,30)]
DEFAULT_HUE = HUES[0]
DEFAULT_SV = [100.0/100.0, 70.0/100.0]
BASELINE_TAG_COLOR = [0, 0, 1]  # HSV
FLOURESCENT_SV = [0.75, 1.0]  # SV for fluorescent colors

def check_plt_backend():
    # reset matplotlib backend for non-interactive plotting
    plt.close('all')
    if plt.get_backend() == 'TkAgg':
        matplotlib.use('Agg')

class LineWidthData(Line2D):
    def __init__(self, *args, **kwargs):
        _lw_data = kwargs.pop('linewidth', 1)
        super().__init__(*args, **kwargs)
        self._lw_data = _lw_data

    def _get_lw(self):
        if self.axes is not None:
            ppd = 72./self.axes.figure.dpi
            trans = self.axes.transData.transform
            return ((trans((1, self._lw_data))-trans((0, 0)))*ppd)[1]
        else:
            return 1

    def _set_lw(self, lw):
        self._lw_data = lw

    _linewidth = property(_get_lw, _set_lw)

def plot_agent(ax, data, color, agent_shape, alpha=1):
    x_center = data['boundary']['location'][0]
    y_center = data['boundary']['location'][1]

    # get color, convert to rgb. Strings are already RGB
    if isinstance(color, str):
        rgb = color
    else:
        rgb = hsv_to_rgb(color)

    if agent_shape == 'rectangle':
        theta = data['boundary']['angle'] / PI * 180 + 90  # rotate 90 degrees to match field
        length = data['boundary']['length']
        width = data['boundary']['width']

        # get bottom left position
        x_offset = (width / 2)
        y_offset = (length / 2)
        theta_rad = math.radians(theta)
        dx = x_offset * math.cos(theta_rad) - y_offset * math.sin(theta_rad)
        dy = x_offset * math.sin(theta_rad) + y_offset * math.cos(theta_rad)

        x = x_center - dx
        y = y_center - dy

        # Create a rectangle
        shape = patches.Rectangle(
            (x, y), width, length,
            angle=theta,
            linewidth=2,
            edgecolor='w',
            alpha=alpha,
            facecolor=rgb,
        )
        ax.add_patch(shape)

    elif agent_shape == 'segment':
        theta = data['boundary']['angle'] / PI * 180 + 90  # rotate 90 degrees to match field
        length = data['boundary']['length']
        width = data['boundary']['width']

        membrane_width = 0.1
        membrane_color = [0, 0, 0]
        radius = width / 2

        # get the two ends
        length_offset = (length / 2) - radius
        theta_rad = math.radians(theta)
        dx = - length_offset * math.sin(theta_rad)
        dy = length_offset * math.cos(theta_rad)

        x1 = x_center - dx
        y1 = y_center - dy
        x2 = x_center + dx
        y2 = y_center + dy

        # segment plot
        membrane = LineWidthData(
            [x1, x2], [y1, y2],
            color=membrane_color,
            linewidth=width,
            alpha=alpha,
            solid_capstyle='round')
        line = LineWidthData(
            [x1, x2], [y1, y2],
            color=rgb,
            alpha=alpha,
            linewidth=width-membrane_width,
            solid_capstyle='round')
        ax.add_line(membrane)
        ax.add_line(line)

    elif agent_shape == 'circle':
        diameter = data['boundary']['diameter']

        # get bottom left position
        radius = (diameter / 2)
        x = x_center - radius
        y = y_center - radius

        # Create a circle
        circle = patches.Circle(
            (x, y), radius, linewidth=1, edgecolor='b', alpha=alpha)
        ax.add_patch(circle)

def plot_agents(
    ax, agents, agent_colors={}, agent_shape='segment', dead_color=None,
    alpha=1
):
    '''
    - ax: the axis for plot
    - agents: a dict with {agent_id: agent_data} and
        agent_data a dict with keys location, angle, length, width
    - agent_colors: dict with {agent_id: hsv color}
    - dead_color: List of 3 floats that define HSV color to use for dead
      cells. Dead cells only get treated differently if this is set.
    - alpha: Alpha value for agents.
    '''
    for agent_id, agent_data in agents.items():
        color = agent_colors.get(agent_id, [DEFAULT_HUE]+DEFAULT_SV)
        if dead_color and 'boundary' in agent_data and 'dead' in agent_data['boundary']:
            if agent_data['boundary']['dead']:
                color = dead_color
        plot_agent(ax, agent_data, color, agent_shape, alpha)
    if len(agents) == 1:
        ax.set_title('1 agent', y=1.1)
    else:
        ax.set_title(f'{len(agents)} agents', y=1.1)

def mutate_color(baseline_hsv):
    mutation = 0.1
    new_hsv = [
        (n + np.random.uniform(-mutation, mutation))
        for n in baseline_hsv]
    # wrap hue around
    new_hsv[0] = new_hsv[0] % 1
    # reflect saturation and value
    if new_hsv[1] > 1:
        new_hsv[1] = 2 - new_hsv[1]
    if new_hsv[2] > 1:
        new_hsv[2] = 2 - new_hsv[2]
    return new_hsv

def color_phylogeny(ancestor_id, phylogeny, baseline_hsv, phylogeny_colors={}):
    """
    get colors for all descendants of the ancestor
    through recursive calls to each generation
    """
    phylogeny_colors.update({ancestor_id: baseline_hsv})
    daughter_ids = phylogeny.get(ancestor_id)
    if daughter_ids:
        for daughter_id in daughter_ids:
            daughter_color = mutate_color(baseline_hsv)
            color_phylogeny(daughter_id, phylogeny, daughter_color)
    return phylogeny_colors

def get_phylogeny_colors_from_names(agent_ids, initial_hues=None):
    '''Get agent colors using phlogeny saved in agent_ids
    This assumes the names use daughter_phylogeny_id() from meta_division
    '''
    initial_hues = initial_hues or {}

    # make phylogeny with {mother_id: [daughter_1_id, daughter_2_id]}
    phylogeny = {agent_id: [] for agent_id in agent_ids}
    for agent1, agent2 in itertools.combinations(agent_ids, 2):
        if agent1 == agent2[0:-1]:
            phylogeny[agent1].append(agent2)
        elif agent2 == agent1[0:-1]:
            phylogeny[agent2].append(agent1)

    # get initial ancestors
    daughters = list(phylogeny.values())
    daughters = set([item for sublist in daughters for item in sublist])
    mothers = set(list(phylogeny.keys()))
    ancestors = list(mothers - daughters)

    # agent colors based on phylogeny
    agent_colors = {agent_id: [] for agent_id in agent_ids}
    for agent_id in ancestors:
        hue = initial_hues.get(agent_id, random.choice(HUES))  # select random initial hue
        initial_color = [hue] + DEFAULT_SV
        agent_colors.update(color_phylogeny(agent_id, phylogeny, initial_color))

    return agent_colors


def plot_snapshots(data, plot_config):
    '''Plot snapshots of the simulation over time

    The snapshots depict the agents and environmental molecule
    concentrations.

    Arguments:
        data (dict): A dictionary with the following keys:

            * **agents** (:py:class:`dict`): A mapping from times to
              dictionaries of agent data at that timepoint. Agent data
              dictionaries should have the same form as the hierarchy
              tree rooted at ``agents``.
            * **fields** (:py:class:`dict`): A mapping from times to
              dictionaries of environmental field data at that
              timepoint.  Field data dictionaries should have the same
              form as the hierarchy tree rooted at ``fields``.
            * **config** (:py:class:`dict`): The environmental
              configuration dictionary  with the following keys:

                * **bounds** (:py:class:`tuple`): The dimensions of the
                  environment.

        plot_config (dict): Accepts the following configuration options.
            Any options with a default is optional.

            * **n_snapshots** (:py:class:`int`): Number of snapshots to
              show per row (i.e. for each molecule). Defaults to 6.
            * **snapshot_times** (:py:class:`Iterable`): Times to plot
              snapshots for. Defaults to None, in which case n_snapshots
              is used.
            * **agent_shape** (:py:class:`str`): the shape of the agents.
              select from **rectangle**, **segment**
            * **phylogeny_names** (:py:class:`bool`): This selects agent
              colors based on phylogenies seved in their names using
              meta_division.py daughter_phylogeny_id()
            * **out_dir** (:py:class:`str`): Output directory, which is
              ``out`` by default.
            * **filename** (:py:class:`str`): Base name of output file.
              ``snapshots`` by default.
            * **skip_fields** (:py:class:`Iterable`): Keys of fields to
              exclude from the plot. This takes priority over
              ``include_fields``.
            * **include_fields** (:py:class:`Iterable`): Keys of fields
              to plot.
            * **field_label_size** (:py:class:`float`): Font size of the
              field label.
            * **dead_color** (:py:class:`list` of 3 :py:class:`float`s):
              Color for dead cells in HSV. Defaults to [0, 0, 0], which
              is black.
            * **default_font_size** (:py:class:`float`): Font size for
              titles and axis labels.
            * **agent_fill_color** (:py:class:`tuple`): HSV Color to use
              for all agent fills.
            * **agent_alpha** (:py:class:`float`): Alpha for agent
              plots.
            * **scale_bar_length** (:py:class:`float`): Length of scale
              bar.  Defaults to 1 (in units of micrometers). If 0, no
              bar plotted.
            * **scale_bar_color** (:py:class:`str`): Color of scale bar
            * **xlim** (:py:class:`tuple` of :py:class:`float`): Tuple
              of lower and upper x-axis limits.
            * **ylim** (:py:class:`tuple` of :py:class:`float`): Tuple
              of lower and upper y-axis limits.
            * **min_color** (any valid matplotlib color): Color for
              minimum field values.
            * **max_color** (any valid matplotlib color): Color for
              maximum field values.
            * **grid_color** (any valid matplotlib color): Whether to
              show. Defaults to an empty string, in which case no
              gridlines are shown.
    '''
    check_plt_backend()

    n_snapshots = plot_config.get('n_snapshots', 6)
    snapshot_times = plot_config.get('snapshot_times')
    out_dir = plot_config.get('out_dir', 'out')
    filename = plot_config.get('filename', 'snapshots')
    agent_shape = plot_config.get('agent_shape', 'segment')
    phylogeny_names = plot_config.get('phylogeny_names', True)
    skip_fields = plot_config.get('skip_fields', [])
    include_fields = plot_config.get('include_fields', None)
    field_label_size = plot_config.get('field_label_size', 20)
    default_font_size = plot_config.get('default_font_size', 36)
    dead_color = plot_config.get('dead_color', [0, 0, 0])
    agent_fill_color = plot_config.get('agent_fill_color', None)
    agent_alpha = plot_config.get('agent_alpha', 1)
    scale_bar_length = plot_config.get('scale_bar_length', 1)
    scale_bar_color = plot_config.get('scale_bar_color', 'black')
    xlim = plot_config.get('xlim')
    ylim = plot_config.get('ylim')
    min_color = plot_config.get('min_color', 'white')
    max_color = plot_config.get('max_color', 'gray')
    grid_color = plot_config.get('grid_color', '')
    agent_colors = plot_config.get('agent_colors', {})
    field_range = plot_config.get('field_range', {})
    begin_gradient = plot_config.get('begin_gradient', 0.25)

    # get data
    agents = data.get('agents', {})
    fields = data.get('fields', {})
    config = data.get('config', {})
    bounds = config.get('bounds', DEFAULT_BOUNDS)
    edge_length_x = bounds[0]
    edge_length_y = bounds[1]

    # time steps that will be used
    if agents and fields:
        assert set(list(agents.keys())) == set(list(fields.keys())), 'agent and field times are different'
        time_vec = list(agents.keys())
    elif agents:
        time_vec = list(agents.keys())
    elif fields:
        time_vec = list(fields.keys())
    else:
        raise Exception('No agents or field data')

    if snapshot_times:
        n_snapshots = len(snapshot_times)
        time_indices = [
            time_vec.index(time) for time in snapshot_times]
    else:
        time_indices = np.round(np.linspace(0, len(time_vec) - 1, n_snapshots)).astype(int)
        snapshot_times = [time_vec[i] for i in time_indices]

    # get fields id and range
    field_ids = []
    if fields:
        if include_fields is None:
            field_ids = set(fields[time_vec[0]].keys())
        else:
            field_ids = set(include_fields)
        field_ids -= set(skip_fields)
        for field_id in field_ids:
            if field_id not in field_range:
                field_min = min([min(min(field_data[field_id])) for t, field_data in fields.items()])
                field_max = max([max(max(field_data[field_id])) for t, field_data in fields.items()])
                field_range[field_id] = [field_min, field_max]

    # get agent ids
    agent_ids = set()
    if agents:
        for time, time_data in agents.items():
            current_agents = list(time_data.keys())
            agent_ids.update(current_agents)
        agent_ids = list(agent_ids)

        # set agent colors
        if agent_colors:
            pass
        elif agent_fill_color:
            agent_colors = {
                agent_id: agent_fill_color for agent_id in agent_ids}
        elif phylogeny_names:
            agent_colors = get_phylogeny_colors_from_names(agent_ids)
        else:
            agent_colors = {}
            for agent_id in agent_ids:
                hue = random.choice(HUES)  # select random initial hue
                color = [hue] + DEFAULT_SV
                agent_colors[agent_id] = color

    # make the figure
    n_rows = max(len(field_ids), 1)
    n_cols = n_snapshots + 1  # one column for the colorbar
    figsize = (12 * n_cols, 12 * n_rows)
    max_dpi = min([2**16 // dim for dim in figsize]) - 1
    fig = plt.figure(figsize=figsize, dpi=min(max_dpi, 100))
    grid = plt.GridSpec(n_rows, n_cols, wspace=0.2, hspace=0.2)
    original_fontsize = plt.rcParams['font.size']
    plt.rcParams.update({'font.size': default_font_size})

    # Add time axis across subplots
    super_spec = matplotlib.gridspec.SubplotSpec(
        grid,
        (n_rows - 1) * n_cols,
        (n_rows - 1) * n_cols + n_snapshots - 1,
    )
    grid_params = grid.get_subplot_params()
    snapshot_times_hrs = tuple(
        time / 60 / 60 for time in snapshot_times)
    if n_snapshots > 1:
        time_per_snapshot = (
            snapshot_times_hrs[-1] - snapshot_times_hrs[0]) / (
            (n_snapshots - 1) * (grid_params.wspace + 1))
    else:
        time_per_snapshot = 1  # Arbitrary
    super_ax = fig.add_subplot(  # type: ignore
        super_spec,
        xticks=snapshot_times_hrs,
        xlim=(
            snapshot_times_hrs[0] - time_per_snapshot / 2,
            snapshot_times_hrs[-1] + time_per_snapshot / 2,
        ),
        yticks=[],
    )
    super_ax.set_xlabel(  # type: ignore
        'Time (hr)', labelpad=50)
    super_ax.xaxis.set_tick_params(width=2, length=8)
    super_ax.xaxis.set_major_formatter(
        matplotlib.ticker.FormatStrFormatter('%.1f'))
    for spine_name in ('top', 'right', 'left'):
        super_ax.spines[spine_name].set_visible(False)
    super_ax.spines['bottom'].set_linewidth(2)

    # Make the colormap
    min_rgb = matplotlib.colors.to_rgb(min_color)
    max_rgb = matplotlib.colors.to_rgb(max_color)
    colors_dict = {
        'red': [
            [0, min_rgb[0], min_rgb[0]],
            [begin_gradient, max_rgb[0], max_rgb[0]],
            [1, max_rgb[0], max_rgb[0]],
        ],
        'green': [
            [0, min_rgb[1], min_rgb[1]],
            [begin_gradient, max_rgb[1], max_rgb[1]],
            [1, max_rgb[1], max_rgb[1]],
        ],
        'blue': [
            [0, min_rgb[2], min_rgb[2]],
            [begin_gradient, max_rgb[2], max_rgb[2]],
            [1, max_rgb[2], max_rgb[2]],
        ],
    }
    cmap = matplotlib.colors.LinearSegmentedColormap(
        'field', segmentdata=colors_dict, N=512)

    stats = {
        'agents': {},
    }
    if field_ids:
        stats['fields'] = {
            field_id: {} for field_id in field_ids
        }
    # plot snapshot data in each subsequent column
    for col_idx, (time_idx, time) in enumerate(zip(time_indices, snapshot_times)):
        stats['agents'][time] = len(agents[time])
        if field_ids:
            for row_idx, field_id in enumerate(field_ids):

                ax = init_axes(
                    fig, edge_length_x, edge_length_y, grid, row_idx,
                    col_idx, time, field_id, field_label_size,
                )
                ax.tick_params(
                    axis='both', which='both', bottom=False, top=False,
                    left=False, right=False,
                )

                # transpose field to align with agents
                field = np.transpose(np.array(fields[time][field_id]))
                vmin, vmax = field_range[field_id]
                q1, q2, q3 = np.percentile(field, [25, 50, 75])
                stats['fields'][field_id][time] = (
                    field.min(), q1, q2, q3, field.max())
                im = plt.imshow(field.tolist(),
                                origin='lower',
                                extent=[0, edge_length_x, 0, edge_length_y],
                                vmin=vmin,
                                vmax=vmax,
                                cmap=cmap)
                ax.set_yticks(
                    np.linspace(0, edge_length_y, field.shape[0] + 1))
                ax.set_xticks(
                    np.linspace(0, edge_length_x, field.shape[1] + 1))
                if agents:
                    agents_now = agents[time]
                    plot_agents(
                        ax, agents_now, agent_colors, agent_shape,
                        dead_color, agent_alpha)
                if grid_color:
                    ax.grid(
                        which='both',
                        color=grid_color,
                        linestyle='-',
                        linewidth=1,
                    )

                if xlim:
                    ax.set_xlim(*xlim)
                if ylim:
                    ax.set_ylim(*ylim)

                # colorbar in new column after final snapshot
                if col_idx == n_snapshots - 1:
                    cbar_col = col_idx + 1
                    cbar_ax = fig.add_subplot(grid[row_idx, cbar_col])
                    if row_idx == 0:
                        cbar_ax.set_title('Concentration (mM)', y=1.08)
                    cbar_ax.axis('off')
                    if vmin == vmax:
                        continue
                    divider = make_axes_locatable(cbar_ax)
                    cax = divider.append_axes("left", size="5%", pad=0.0)
                    fig.colorbar(im, cax=cax, format='%.3f')
                    cbar_ax.axis('off')
                # Scale bar in first snapshot of each row
                if col_idx == 0 and scale_bar_length:
                    scale_bar = anchored_artists.AnchoredSizeBar(
                        ax.transData, scale_bar_length,
                        f'${scale_bar_length} \\mu m$', 'lower left',
                        color=scale_bar_color,
                        frameon=False,
                        sep = scale_bar_length,
                        size_vertical = scale_bar_length / 20,
                    )
                    ax.add_artist(scale_bar)
        else:
            row_idx = 0
            ax = init_axes(
                fig, bounds[0], bounds[1], grid, row_idx, col_idx,
                time, ""
            )
            ax.tick_params(
                axis='both', which='both', bottom=False, top=False,
                left=False, right=False,
            )

            if agents:
                agents_now = agents[time]
                plot_agents(
                    ax, agents_now, agent_colors, agent_shape,
                    dead_color, agent_alpha)
            if xlim:
                ax.set_xlim(*xlim)
            if ylim:
                ax.set_ylim(*ylim)
            # Scale bar in first snapshot of each row
            if col_idx == 0 and scale_bar_length:
                scale_bar = anchored_artists.AnchoredSizeBar(
                    ax.transData, scale_bar_length,
                    f'${scale_bar_length} \\mu m$', 'lower left',
                    color=scale_bar_color,
                    frameon=False,
                    sep = scale_bar_length,
                    size_vertical = scale_bar_length / 20,
                )
                ax.add_artist(scale_bar)

    fig_path = os.path.join(out_dir, filename)
    fig.subplots_adjust(wspace=0.7, hspace=0.1)
    fig.savefig(fig_path, bbox_inches='tight')
    plt.close(fig)
    plt.rcParams.update({'font.size': original_fontsize})
    return stats


def get_tag_ranges(agents, tagged_molecules, convert_to_concs):
    tag_ranges = {}
    for time, time_data in agents.items():
        for agent_id, agent_data in time_data.items():
            volume = agent_data.get('boundary', {}).get('volume', 0)
            for tag_id in tagged_molecules:
                level = get_value_from_path(agent_data, tag_id)
                if convert_to_concs:
                    level = level / volume if volume else 0
                if tag_id in tag_ranges:
                    tag_ranges[tag_id] = [
                        min(tag_ranges[tag_id][0], level),
                        max(tag_ranges[tag_id][1], level)]
                else:
                    # add new tag
                    tag_ranges[tag_id] = [level, level]
    return tag_ranges


def plot_tags(data, plot_config):
    '''Plot snapshots of the simulation over time

    The snapshots depict the agents and the levels of tagged molecules
    in each agent by agent color intensity.

    Arguments:
        data (dict): A dictionary with the following keys:

            * **agents** (:py:class:`dict`): A mapping from times to
              dictionaries of agent data at that timepoint. Agent data
              dictionaries should have the same form as the hierarchy
              tree rooted at ``agents``.
            * **config** (:py:class:`dict`): The environmental
              configuration dictionary  with the following keys:

                * **bounds** (:py:class:`tuple`): The dimensions of the
                  environment.

        plot_config (dict): Accepts the following configuration options.
            Any options with a default is optional.

            * **n_snapshots** (:py:class:`int`): Number of snapshots to
              show per row (i.e. for each molecule). Defaults to 6.
            * **out_dir** (:py:class:`str`): Output directory, which is
              ``out`` by default.
            * **filename** (:py:class:`str`): Base name of output file.
              ``tags`` by default.
            * **tagged_molecules** (:py:class:`typing.Iterable`): The
              tagged molecules whose concentrations will be indicated by
              agent color. Each molecule should be specified as a
              :py:class:`tuple` of the path in the agent compartment
              to where the molecule's count can be found, with the last
              value being the molecule's count variable.
            * **convert_to_concs** (:py:class:`bool`): if True, convert counts
              to concentrations.
            * **background_color** (:py:class:`str`): use matplotlib colors,
              ``black`` by default
            * **tag_label_size** (:py:class:`float`): The font size for
              the tag name label
            * **default_font_size** (:py:class:`float`): Font size for
              titles and axis labels.
            * **tag_colors** (:py:class:`dict`): Map from tag ID in
              tagged_molecules to a tuple (min_color, max_color).
            * **scale_bar_length** (:py:class:`float`): Length of scale
              bar.  Defaults to 1 (in units of micrometers). If 0, no
              bar plotted.
            * **scale_bar_color** (:py:class:`str`): Color of scale bar
            * **xlim** (:py:class:`tuple` of :py:class:`float`): Tuple
              of lower and upper x-axis limits.
            * **ylim** (:py:class:`tuple` of :py:class:`float`): Tuple
              of lower and upper y-axis limits.
    '''
    check_plt_backend()

    n_snapshots = plot_config.get('n_snapshots', 6)
    out_dir = plot_config.get('out_dir', 'out')
    filename = plot_config.get('filename', 'tags')
    agent_shape = plot_config.get('agent_shape', 'segment')
    background_color = plot_config.get('background_color', 'black')
    tagged_molecules = plot_config['tagged_molecules']
    tag_path_name_map = plot_config.get('tag_path_name_map', {})
    tag_label_size = plot_config.get('tag_label_size', 20)
    default_font_size = plot_config.get('default_font_size', 36)
    convert_to_concs = plot_config.get('convert_to_concs', True)
    tag_colors = plot_config.get('tag_colors', {})
    scale_bar_length = plot_config.get('scale_bar_length', 1)
    scale_bar_color = plot_config.get('scale_bar_color', 'white')
    xlim = plot_config.get('xlim')
    ylim = plot_config.get('ylim')
    tag_ranges = plot_config.get('tag_ranges', {})
    snapshot_times = plot_config.get('snapshot_times', [])

    if tagged_molecules == []:
        raise ValueError('At least one molecule must be tagged.')

    # get data
    agents = data['agents']
    config = data.get('config', {})
    bounds = config['bounds']
    edge_length_x, edge_length_y = bounds

    # time steps that will be used
    time_vec = list(agents.keys())
    if snapshot_times:
        n_snapshots = len(snapshot_times)
        time_indices = [
            time_vec.index(time) for time in snapshot_times]
    else:
        time_indices = np.round(np.linspace(0, len(time_vec) - 1, n_snapshots)).astype(int)
        snapshot_times = [time_vec[i] for i in time_indices]

    # get tag ids and range
    if not tag_ranges:
        tag_ranges = get_tag_ranges(
            agents, tagged_molecules, convert_to_concs)

    # make the figure
    n_rows = len(tagged_molecules)
    n_cols = n_snapshots + 1  # one column for the colorbar
    figsize = (12 * n_cols, 12 * n_rows)
    max_dpi = min([2**16 // dim for dim in figsize]) - 1
    fig = plt.figure(figsize=figsize, dpi=min(max_dpi, 100))
    grid = plt.GridSpec(n_rows, n_cols, wspace=0.2, hspace=0.2)
    original_fontsize = plt.rcParams['font.size']
    plt.rcParams.update({'font.size': default_font_size})

    # Add time axis across subplots
    super_spec = matplotlib.gridspec.SubplotSpec(
        grid,
        (n_rows - 1) * n_cols,
        (n_rows - 1) * n_cols + n_snapshots - 1,
    )
    grid_params = grid.get_subplot_params()
    snapshot_times_hrs = tuple(time / 60 / 60 for time in snapshot_times)
    if n_snapshots > 1:
        time_per_snapshot = (
            snapshot_times_hrs[-1] - snapshot_times_hrs[0]) / (
            (n_snapshots - 1) * (grid_params.wspace + 1))
    else:
        time_per_snapshot = 1  # Arbitrary
    super_ax = fig.add_subplot(  # type: ignore
        super_spec,
        xticks=snapshot_times_hrs,
        xlim=(
            snapshot_times_hrs[0] - time_per_snapshot / 2,
            snapshot_times_hrs[-1] + time_per_snapshot / 2,
        ),
        yticks=[],
    )
    super_ax.set_xlabel(  # type: ignore
        'Time (hr)', labelpad=50)
    super_ax.xaxis.set_tick_params(width=2, length=8)
    super_ax.xaxis.set_major_formatter(
        matplotlib.ticker.FormatStrFormatter('%.1f'))
    for spine_name in ('top', 'right', 'left'):
        super_ax.spines[spine_name].set_visible(False)
    super_ax.spines['bottom'].set_linewidth(2)

    # plot tags
    for row_idx, tag_id in enumerate(tag_ranges.keys()):
        tag_name = tag_path_name_map.get(tag_id, tag_id)
        min_tag, max_tag = tag_ranges[tag_id]
        min_color, max_color = tag_colors[tag_id]
        min_rgb = matplotlib.colors.to_rgb(min_color)
        max_rgb = matplotlib.colors.to_rgb(max_color)
        colors_dict = {
            'red': [
                [0, min_rgb[0], min_rgb[0]],
                [1, max_rgb[0], max_rgb[0]],
            ],
            'green': [
                [0, min_rgb[1], min_rgb[1]],
                [1, max_rgb[1], max_rgb[1]],
            ],
            'blue': [
                [0, min_rgb[2], min_rgb[2]],
                [1, max_rgb[2], max_rgb[2]],
            ],
        }
        cmap = matplotlib.colors.LinearSegmentedColormap(
            tag_id, segmentdata=colors_dict, N=512)

        norm = matplotlib.colors.Normalize(min_tag, max_tag)

        for col_idx, (_, time) in enumerate(
            zip(time_indices, snapshot_times)
        ):
            ax = init_axes(
                fig, edge_length_x, edge_length_y, grid,
                row_idx, col_idx, time, tag_name, tag_label_size,
            )
            ax.tick_params(
                axis='both', which='both', bottom=False, top=False,
                left=False, right=False,
            )
            ax.set_facecolor(background_color)

            agent_tag_colors = {}
            for agent_id, agent_data in agents[time].items():
                # get current tag concentration, and determine color
                level = get_value_from_path(agent_data, tag_id)
                if convert_to_concs:
                    volume = get_in(
                        agent_data, ('boundary', 'volume'), 0)
                    level = level / volume if volume else 0
                intensity = norm(level)
                agent_rgb = cmap(intensity)[:3]
                agent_hsv = matplotlib.colors.rgb_to_hsv(
                    agent_rgb)
                agent_tag_colors[agent_id] = agent_hsv

            plot_agents(ax, agents[time], agent_tag_colors, agent_shape)

            if xlim:
                ax.set_xlim(*xlim)
            if ylim:
                ax.set_ylim(*ylim)

            # colorbar in new column after final snapshot
            if col_idx == n_snapshots - 1:
                cbar_col = col_idx + 1
                cbar_ax = fig.add_subplot(grid[row_idx, cbar_col])
                if row_idx == 0:
                    if convert_to_concs:
                        cbar_ax.set_title('Concentration (counts/fL)', y=1.08)
                cbar_ax.axis('off')
                if min_tag == max_tag:
                    continue
                divider = make_axes_locatable(cbar_ax)
                cax = divider.append_axes("left", size="5%", pad=0.0)
                mappable = matplotlib.cm.ScalarMappable(norm, cmap)
                fig.colorbar(mappable, cax=cax, format='%.0f')

            # Scale bar in first snapshot of each row
            if col_idx == 0 and scale_bar_length:
                scale_bar = anchored_artists.AnchoredSizeBar(
                    ax.transData, scale_bar_length,
                    f'${scale_bar_length} \\mu m$', 'lower left',
                    color=scale_bar_color,
                    frameon=False,
                    sep = scale_bar_length,
                    size_vertical = scale_bar_length / 20,
                )
                ax.add_artist(scale_bar)

    fig_path = os.path.join(out_dir, filename)
    fig.subplots_adjust(wspace=0.7, hspace=0.1)
    fig.savefig(fig_path, bbox_inches='tight')
    plt.close(fig)
    plt.rcParams.update({'font.size': original_fontsize})




def init_axes(
    fig, edge_length_x, edge_length_y, grid, row_idx, col_idx, time,
    molecule, ylabel_size=20
):
    ax = fig.add_subplot(grid[row_idx, col_idx])
    if row_idx == 0:
        plot_title = 'time: {:.0f} s'.format(float(time))
        plt.title(plot_title, y=1.08)
    if col_idx == 0:
        ax.set_ylabel(
            molecule, fontsize=ylabel_size, rotation='horizontal',
            horizontalalignment='right', labelpad=50,
        )
    ax.set(xlim=[0, edge_length_x], ylim=[0, edge_length_y], aspect=1)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    return ax
