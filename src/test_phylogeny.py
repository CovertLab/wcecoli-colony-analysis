from typing import Tuple, cast

from ete3 import TreeNode

from src.phylogeny import make_ete_trees


class TestMakeEteTree:

    @staticmethod
    def _traverse(tree: TreeNode) -> Tuple[str]:
        traversal = tuple(
            node.name
            for node in tree.traverse('preorder')
        )
        return cast(Tuple[str], traversal)

    def test_simple(self) -> None:
        # Tree:
        #                   /-agent00
        #         /-agent0-|
        # agent--|          \-agent01
        #        |
        #         \-agent1
        agent_ids = (
            'agent', 'agent0', 'agent1', 'agent00', 'agent01')
        trees = make_ete_trees(agent_ids)
        assert len(trees) == 1
        traversal = self._traverse(trees[0])
        expected_traversal = (
            'agent', 'agent0', 'agent00', 'agent01', 'agent1')
        assert traversal == expected_traversal

    @staticmethod
    def test_empty() -> None:
        trees = make_ete_trees([])
        assert trees == []

    def test_single(self) -> None:
        trees = make_ete_trees(('agent',))
        assert len(trees) == 1
        traversal = self._traverse(trees[0])
        expected_traversal = ('agent',)
        assert traversal == expected_traversal
