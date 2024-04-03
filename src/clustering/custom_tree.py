from ete3 import Tree


def read_tree(tree_path: str):
    tree = Tree(tree_path, format=5)

    # Access the tree properties
    print("Root node:", tree.get_tree_root())
    print("Number of leaves:", tree.get_leaf_count())

    # Traverse the tree
    for node in tree.traverse():
        print("Node name:", node.name)

    # Access parent-child relationships
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            print("Leaf:", node.name)
        else:
            children = [child.name for child in node.children]
            print("Internal node:", node.name, "Children:", children)
