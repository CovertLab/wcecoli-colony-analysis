import math
import os
import shutil

import cv2  # type: ignore
from vivarium_cell.analysis.analyze import Analyzer

from src.plot_snapshots import (
    plot_snapshots,
    plot_tags,
    get_phylogeny_colors_from_names,
    get_tag_ranges,
)
from src.types import RawData


# 6.4 hours of simulation time in 10 seconds
SPEED_SCALE = 10 / (6.4 * 60 * 60)
FPS = 15  # Frames-per-second of final video


def plot_single_snapshot(
        data, time, fields, out_dir, filename, agent_fill_color, xlim,
        ylim, agent_alpha, agent_colors, field_range, scalebar_color):
    plot_config = {
        'out_dir': out_dir,
        'filename': filename,
        'include_fields': fields,
        'field_label_size': 36,
        'default_font_size': 36,
        'agent_fill_color': agent_fill_color,
        'dead_color': (0, 0, 0.79),  # gray in HSV
        'agent_alpha': agent_alpha,
        'n_snapshots': 1,
        'snapshot_times': [time],
        'scale_bar_length': 10,
        'scale_bar_color': scalebar_color,
        'xlim': xlim,
        'ylim': ylim,
        'min_color': '#FFFFFF',
        'max_color': '#000000',
        'grid_color': 'white' if fields else '',
        'agent_colors': agent_colors,
        'field_range': field_range,
        'begin_gradient': 1,
    }
    plot_snapshots(data, plot_config)


def plot_single_tags_plot(
        data, time, tags, tag_path_map, tag_ranges, out_dir, filename,
        xlim, ylim):
    plot_config = {
        'out_dir': out_dir,
        'tagged_molecules': tags,
        'background_color': 'white',
        'filename': filename,
        'tag_path_name_map': tag_path_map,
        'tag_label_size': 36,
        'default_font_size': 36,
        'n_snapshots': 1,
        'tag_colors': {
            tag: ('white', '#0000ff')
            for tag in tags
        },
        'scale_bar_length': 10,
        'scale_bar_color': 'black',
        'xlim': xlim,
        'ylim': ylim,
        'snapshot_times': [time],
        'tag_ranges': tag_ranges,
    }
    plot_tags(data, plot_config)


def make_tags_video(
        data,
        environment_config,
        tags,
        tag_path_map,
        out_dir='out',
        filename='snapshot_vid',
        xlim=None,
        ylim=None,
        ):
    # make images directory, remove if existing
    out_file = os.path.join(out_dir, f'{filename}.mp4')
    images_dir = os.path.join(out_dir, '_images')
    os.makedirs(images_dir)

    all_times = list(data.keys())
    max_time = max(all_times)
    num_frames = len(all_times)
    desired_runtime = max_time * SPEED_SCALE
    desired_num_frames = desired_runtime * FPS
    available_frames_per_desired = math.ceil(
        num_frames / desired_num_frames)
    time_vec = all_times[::available_frames_per_desired]
    true_fps = len(time_vec) / desired_runtime

    tags_data = Analyzer.format_data_for_tags(data, environment_config)
    tag_ranges = get_tag_ranges(tags_data['agents'], tags, True)

    # make the individual snapshot figures
    img_paths = []
    for i, time in enumerate(time_vec):
        img_file = f'img_{i}.png'
        plot_single_tags_plot(
            tags_data, time, tags, tag_path_map, tag_ranges, images_dir,
            img_file, xlim, ylim,
        )
        fig_path = os.path.join(images_dir, img_file)
        img_paths.append(fig_path)

    # make the video
    img_array = []
    for img_file in img_paths:
        img = cv2.imread(img_file)
        height, width, _ = img.shape
        size = (width, height)
        img_array.append(img)

    out = cv2.VideoWriter(
        out_file, cv2.VideoWriter_fourcc(*'mp4v'), true_fps, size)

    for img_path in img_array:
        out.write(img_path)
    out.release()

    # delete image folder
    shutil.rmtree(images_dir)


def make_snapshots_video(
        data,
        environment_config,
        fields,
        out_dir='out',
        filename='snapshot_vid',
        xlim=None,
        ylim=None,
        agent_alpha=1,
        agent_fill_color=None,
        field_added=0,
        ):
    # make images directory, remove if existing
    out_file = os.path.join(out_dir, f'{filename}.mp4')
    images_dir = os.path.join(out_dir, '_images')
    os.makedirs(images_dir)

    all_times = list(data.keys())
    max_time = max(all_times)
    num_frames = len(all_times)
    desired_runtime = max_time * SPEED_SCALE
    desired_num_frames = desired_runtime * FPS
    available_frames_per_desired = math.ceil(
        num_frames / desired_num_frames)
    time_vec = all_times[::available_frames_per_desired]
    true_fps = len(time_vec) / desired_runtime

    snapshots_data = Analyzer.format_data_for_snapshots(
        data, environment_config)

    agent_ids = set()
    for time_data in snapshots_data['agents'].values():
        agent_ids.update(time_data.keys())

    # get fields and agent colors
    field_range = {}
    fields_data = snapshots_data['fields']
    for field in fields:
        field_min = min([
            min(min(field_data[field]))
            for t, field_data in fields_data.items()
        ])
        field_max = max([
            max(max(field_data[field]))
            for t, field_data in fields_data.items()
        ])
        field_range[field] = [field_min, field_max]
    agent_colors = get_phylogeny_colors_from_names(agent_ids)

    if not fields:
        data = RawData({
            key: val
            for key, val in data.items() if key != 'fields'
        })

    # make the individual snapshot figures
    img_paths = []
    for i, time in enumerate(time_vec):
        scalebar_color = (
            'black' if time < max_time * field_added else 'white')
        img_file = f'img_{i}.png'
        plot_single_snapshot(
            snapshots_data, time, fields, images_dir, img_file,
            agent_fill_color, xlim, ylim, agent_alpha, agent_colors,
            field_range, scalebar_color,
        )
        fig_path = os.path.join(images_dir, img_file)
        img_paths.append(fig_path)

    # make the video
    img_array = []
    for img_file in img_paths:
        img = cv2.imread(img_file)
        height, width, _ = img.shape
        size = (width, height)
        img_array.append(img)

    out = cv2.VideoWriter(
        out_file, cv2.VideoWriter_fourcc(*'mp4v'), true_fps, size)

    for img_path in img_array:
        out.write(img_path)
    out.release()

    # delete image folder
    shutil.rmtree(images_dir)
