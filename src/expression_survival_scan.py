import json
from matplotlib import pyplot as plt

from vivarium.core.experiment import get_in


LIVE_COLOR = 'red'
DEAD_COLOR = 'black'
ALPHA = 0.5


def plot_expression_survival_scan(raw_datasets, agent_name):
    death_path = ('agents', agent_name, 'boundary', 'dead')
    points_die = []  # (pump, beta lactamase)
    points_live = []  # (pump, beta lactamase)
    for pump, beta_lactamase, raw_data in raw_datasets:
        max_time = max(float(key) for key in raw_data.keys())
        last_timepoint = raw_data[str(max_time)]
        dead = get_in(last_timepoint, death_path)
        if dead:
            points_die.append((pump, beta_lactamase))
        else:
            points_live.append((pump, beta_lactamase))

    fig, ax = plt.subplots()
    if points_die:
        ax.scatter(
            *zip(*points_die), label='Die', color=DEAD_COLOR,
            alpha=ALPHA)
    if points_live:
        ax.scatter(
            *zip(*points_live), label='Survive', color=LIVE_COLOR,
            alpha=ALPHA)
    ax.legend()
    ax.set_xlabel('[AcrAB-TolC] (mM)')
    ax.set_ylabel('[AmpC] (mM)')
    for spine_name in ('top', 'right'):
        ax.spines[spine_name].set_visible(False)
    fig.tight_layout()
    return fig


def load_scan_data(path):
    with open(path, 'r') as f:
        data_dict = json.load(f)
    return data_dict['data'], data_dict['parameters']
