import json

from matplotlib import pyplot as plt
import numpy as np
from sklearn.linear_model import LogisticRegression

from vivarium.core.experiment import get_in


LIVE_COLOR = 'green'
DEAD_COLOR = 'black'
ALPHA = 0.5


def plot_expression_survival_scan(
        raw_datasets, agent_name, xlabel, ylabel, numeric_x, numeric_y,
        numeric_error, scaling=1):
    death_path = ('agents', agent_name, 'boundary', 'dead')
    points_die = []  # (pump, beta lactamase)
    points_live = []  # (pump, beta lactamase)
    for pump, beta_lactamase, raw_data in raw_datasets:
        max_time = max(float(key) for key in raw_data.keys())
        last_timepoint = raw_data[str(max_time)]
        dead = get_in(last_timepoint, death_path)
        if dead:
            points_die.append(
                (pump * scaling, beta_lactamase * scaling))
        else:
            points_live.append(
                (pump * scaling, beta_lactamase * scaling))

    fig, ax = plt.subplots()
    if points_die:
        ax.scatter(
            *zip(*points_die), label='Die', color=DEAD_COLOR,
            alpha=ALPHA)
    if points_live:
        ax.scatter(
            *zip(*points_live), label='Survive', color=LIVE_COLOR,
            alpha=ALPHA)

    model = LogisticRegression()
    features = np.array(points_die + points_live)
    # Logistic regression doesn't do well with very small numbers
    scaled = features / features.max()
    labels = [0] * len(points_die) + [1] * len(points_live)
    model.fit(scaled, labels)
    theta_0 = model.intercept_[0]
    theta_1, theta_2 = model.coef_.tolist()[0]
    x, _ = zip(*features)
    boundary_x = np.linspace(min(x), max(x), 10)
    # Decision boundary comes from solving theta^T X = 0
    m = -theta_1 / theta_2
    b = -theta_0 / theta_2 * features.max()
    boundary_y = m * boundary_x + b
    ax.plot(
        boundary_x, boundary_y, c='blue',
        label='y = {}x + {}'.format(m, b))

    # Plot numeric solution
    numeric_x_scaled = np.array(numeric_x) * scaling
    numeric_y_scaled = np.array(numeric_y) * scaling
    numeric_error_scaled = numeric_error * scaling
    ax.plot(
        numeric_x_scaled, numeric_y_scaled, c='orange',
        label='Numeric Solution for Boundary')
    ax.fill_between(
        numeric_x_scaled,
        numeric_y_scaled - numeric_error_scaled / 2,
        numeric_y_scaled + numeric_error_scaled / 2,
        color='orange', alpha=0.2
    )

    ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    for spine_name in ('top', 'right'):
        ax.spines[spine_name].set_visible(False)
    fig.tight_layout()
    return fig


def load_scan_data(path):
    with open(path, 'r') as f:
        data_dict = json.load(f)
    return data_dict['data'], data_dict['parameters']
