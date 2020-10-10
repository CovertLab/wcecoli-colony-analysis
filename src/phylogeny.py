'''Tools for working with agent phylogenies.'''

import os
from typing import Dict, List, Set, Iterable

from ete3 import TreeNode, TreeStyle, NodeStyle
from vivarium.core.experiment import get_in

from src.types import RawData
from src.constants import AGENTS_PATH


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


def plot_phylogeny(data: RawData, out: str = 'phylogeny.pdf') -> None:
    '''Plot phylogenetic tree from an experiment.

    Args:
        data: The simulation data.
        out: Path to the output file. File type will be inferred from
            the file name.
    '''
    agent_ids: Set[str] = set()
    for _, time_data in data.items():
        agents_data = get_in(time_data, AGENTS_PATH)
        assert agents_data is not None
        agent_ids |= set(agents_data.keys())
    trees = make_ete_trees(agent_ids)
    assert len(trees) == 1
    tree = trees[0]
    tstyle = TreeStyle()
    tstyle.show_scale = False
    tstyle.show_leaf_name = False
    tstyle.scale = 10
    nstyle=NodeStyle()
    nstyle['size'] = 5
    nstyle['vt_line_width'] = 3
    nstyle['hz_line_width'] = 3
    for node in tree.traverse():
        node.set_style(nstyle)
    tree.render(out, tree_style=tstyle, w=400)
