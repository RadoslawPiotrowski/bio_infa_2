import dendropy
from Bio import Phylo
from dendropy import Tree
from dendropy.calculate import treecompare

PATH_TO_TREES = 'nwk_trees'
PATH_EXAMPLE = 'nwk_trees/example.nwk'
PATH_BIF = 'nwk_trees/bifurcating.nwk'
PATH_PRE = 'nwk_trees/preterminal.nwk'


def tree_info(tree) -> str:

    info = f"NUMBER OF NODES: { len(tree.depths()) } \n" \
        f"ALL TERMINALS/LEAFS/NODES: { tree.get_terminals() } \n" \
        f"ALL DEPTHS FROM ROOT: { tree.depths() } \n" \
        f"ALL INTERNAL TERMINALS/LEAFS/NODES: { tree.get_nonterminals() } \n" \
        f"NUMBER OF TERMINALS/LEAFS/NODES: { tree.count_terminals() } \n" \
        f"IS TREE BIFURCATING/BINARY NODES HAS 0 or 2 CHILDS NOT INCLUDE ROOT: { tree.is_bifurcating() } \n" \
        f"IS TREE PRETERMINAL HAS ONLY ONE LEVEL DEPTH: { tree.is_preterminal() } \n" \
        f"TOTAL TREE BRANCH LENGTH: { tree.total_branch_length() } \n"
    return info


def main():
    # getting the tree
    tree_gen = Phylo.parse(PATH_EXAMPLE, 'newick')
    tree_object = next(tree_gen)

    # the tree basic information
    print(tree_info(tree_object))

    # drawing the tree
    Phylo.draw(tree_object)

    # distance comparing
    tns = dendropy.TaxonNamespace()
    tre_one = Tree.get_from_path(PATH_EXAMPLE, 'newick', taxon_namespace=tns)
    tre_two = Tree.get_from_path(PATH_BIF, 'newick', taxon_namespace=tns)

    euclidean_distance = treecompare.euclidean_distance(tre_one, tre_two)
    robinson_distance = treecompare.robinson_foulds_distance(tre_one, tre_two)
    print("Robinson Foulds distance: ", robinson_distance)
    print("Euclidean distance: ", euclidean_distance)

    # common ancestors
    common_ancestor_tree = tree_object.common_ancestor({"name": "C"},
                                    {"name": "D"})
    common_ancestor_tree.color = "blue"
    print("COMMON ANCESTOR: ", common_ancestor_tree)
    Phylo.draw(common_ancestor_tree)


if __name__ == "__main__":
    main()
