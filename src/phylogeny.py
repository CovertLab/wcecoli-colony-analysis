'''Tools for working with agent phylogenies.'''

import os
from typing import Dict, List, Set, Iterable, Tuple

from ete3 import TreeNode, TreeStyle, NodeStyle, CircleFace, TextFace
import pandas as pd
from vivarium.core.experiment import get_in

from src.types import RawData
from src.constants import AGENTS_PATH


PATH_TO_DEAD = ('boundary', 'dead')


def make_ete_trees(agent_ids: Iterable[str]) -> List[TreeNode]:
    '''Construct an ETE Toolkit Tree from a sequence of agent IDs

    Agent IDs must be constructed such that for any agent with ID
    :math:`p` with a parent with ID :math:`p`, :math:`p == c[:-1]`. This
    function should be able to handle multiple phylogenies among the
    agents, but this behavior is not guaranteed, tested, nor supported.

    Args:
        agent_ids: Sequence of agent IDs to build a tree from.

    Returns:
        A list of the roots of the created trees.
    '''
    stem = os.path.commonprefix(list(agent_ids))
    id_node_map: Dict[str, TreeNode] = dict()
    sorted_agents = sorted(agent_ids)
    roots: List[TreeNode] = []
    for agent_id in sorted_agents:
        phylogeny_id = agent_id[len(stem):]
        try:
            if phylogeny_id:
                int(phylogeny_id)
        except ValueError as e:
            raise ValueError(
                'String in ID {} after stem {} is non-numeric'.format(
                    agent_id, stem)
            ) from e
        parent_phylo_id = phylogeny_id[:-1]
        if parent_phylo_id in id_node_map:
            parent = id_node_map[parent_phylo_id]
            child = parent.add_child(name=agent_id)
        else:
            child = TreeNode(name=agent_id)
            roots.append(child)
        id_node_map[phylogeny_id] = child
    return roots


def plot_phylogeny(
        data: RawData, out: str = 'phylogeny.pdf',
        live_color: str = 'green', dead_color: str = 'black',
        ignore_color: str = 'lightgray',
        time_range: Tuple[float, float] = (0, 1)
        ) -> Tuple[TreeNode, pd.DataFrame]:
    '''Plot phylogenetic tree from an experiment.

    Args:
        data: The simulation data.
        out: Path to the output file. File type will be inferred from
            the file name.
        live_color: Color for nodes representing cells that survive
            until division.
        dead_color: Color for nodes representing cells that die.
        ignore_color: Color for nodes outside the time range considered.
        time_range: Tuple specifying the range of times to consider.
            Range values specified as fractions of the final
            timepointpoint.
    '''
    agent_ids: Set[str] = set()
    dead_ids: Set[str] = set()
    in_time_range_ids: Set[str] = set()
    end_time = max(data.keys())
    for time, time_data in data.items():
        agents_data = get_in(time_data, AGENTS_PATH)
        assert agents_data is not None
        agent_ids |= set(agents_data.keys())

        if time_range[0] * end_time <= time <= time_range[1] * end_time:
            in_time_range_ids |= set(agents_data.keys())
            for agent_id, agent_data in agents_data.items():
                if get_in(agent_data, PATH_TO_DEAD, False):
                    dead_ids.add(agent_id)

    trees = make_ete_trees(agent_ids)
    assert len(trees) == 1
    tree = trees[0]

    # Set style for overall figure
    tstyle = TreeStyle()
    tstyle.show_scale = False
    tstyle.show_leaf_name = False
    tstyle.scale = None
    tstyle.optimal_scale_level = 'full'  # Avoid artificial branches
    tstyle.mode = 'c'
    legend = {
        'Die': dead_color,
        'Survive': live_color,
        'Divided Before Antibiotics Appeared': ignore_color,
    }
    for label, color in legend.items():
        tstyle.legend.add_face(CircleFace(5, color), column=0)
        tstyle.legend.add_face(TextFace(label), column=1)

    # Set styles for each node
    for node in tree.traverse():
        nstyle=NodeStyle()
        nstyle['size'] = 5
        nstyle['vt_line_width'] = 1
        nstyle['hz_line_width'] = 1
        if node.name in in_time_range_ids:
            if node.name in dead_ids:
                nstyle['fgcolor'] = dead_color
            else:
                nstyle['fgcolor'] = live_color
        else:
            nstyle['fgcolor'] = ignore_color
        node.set_style(nstyle)
    tree.render(out, tree_style=tstyle, w=400)
    survive_col = []
    agents_col = []
    for agent in in_time_range_ids:
        agents_col.append(agent)
        survive_col.append(0 if agent in dead_ids else 1)
    df = pd.DataFrame({'agents': agents_col, 'survival': survive_col})
    return tree, df
