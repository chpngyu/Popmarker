# ranking genes
#
# Chun-Ping Yu  chpngyu@gmail.com
# July 18, 2016
#

from Bio import Phylo
from matplotlib import pyplot as plt
import matplotlib
from itertools import cycle

class Tree:
    def __init__(self, file_path, tree_format):
        self.tree = Phylo.read(file_path, tree_format)
        self.colors = cycle(['#DF3A01', '#3ADF00', '#0101DF', '#DF01A5', '#DBA901',
                       '#01A9DB', '#21610B', '#610B4B', '#8A4B08', '#0B4C5F'])
        clade_uid = 0
        # assigning clades which have no names
        for clade in self.tree.find_clades():
            if clade.name == None:
                clade.name = 'CID%d' % clade_uid
                clade_uid += 1
        print('There are %d nodes in %s' % (self.tree.count_terminals(), file_path))

    def grouping(self, max_depth, min_size):
        _clade_depth = self.tree.depths(unit_branch_lengths=True)
        # self.depth_to_node stores the nodes which have same depths
        # self.depth_to_node = {3: set([node1, node2]), 4: set([node3, node4, node5]), ...}
        self.depth_to_node = {}
        for clade in _clade_depth:
            if not clade.is_terminal():
                node, depth = clade.name, _clade_depth[clade]
                if depth not in self.depth_to_node:
                    self.depth_to_node[depth] = set([node,])
                else:
                    self.depth_to_node[depth].add(node)
        used_clade = []
        self.grouped_nodes = []
        for depth in range(max_depth, 0, -1):
            grouped_nodes, more_used_clade = self.get_unused_clade(self.depth_to_node[depth], used_clade, min_size)
            if grouped_nodes:
                self.grouped_nodes += grouped_nodes
                used_clade += more_used_clade
        if (depth == 0) and (len(self.grouped_sets) == 0):
            raise RuntimeError('No more grouping clades can be found such that the minimal sizes >= %d' % min_size)            

    def get_unused_clade(self, targets, used_clades, min_size):
        _grouped_sets = []
        _used_clades = list(used_clades)
        for target in targets:
            not_used = True
            clade = next(self.tree.find_clades({'name': target}))
            for uclade in _used_clades:
                common_ancestor = self.tree.common_ancestor({'name': clade.name}, {'name': uclade.name})
                if clade.name == common_ancestor.name:
                    not_used = False
                    break
            if not_used and (clade.count_terminals() >= min_size):
                clade.color = next(self.colors)
                _grouped_sets.append([x.name for x in clade.get_terminals()])
                _used_clades.append(clade)
        return _grouped_sets, _used_clades

    def distance(self, targets):
        _dists = []
        for ix, t1 in enumerate(targets):
            for t2 in targets[ix+1:]:
                _dists.append( self.tree.distance({'name': t1}, {'name': t2}) )
        return _dists

    def save(self, file_path, font_size=6):
        from copy import deepcopy
        _tree = deepcopy(self.tree)
        # remove the clade IDs
        for clade in _tree.find_clades():
            if (not clade.is_terminal()) and (clade.name.find("CID") == 0):
                clade.name = None
        matplotlib.rc('font', **{'size': font_size})
        _tree.ladderize()
        Phylo.draw(_tree, do_show=False)
        plt.savefig('tree.pdf')


if __name__ == '__main__':
    import argparse
    detailed_info = 'Output a (trimmed) tree to a figure'
    parser = argparse.ArgumentParser(description='Separating a tree into smaller trees',
                                     epilog=detailed_info, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--tree', type=str, required=True,
                        help='phylogenetic tree')
    parser.add_argument('--maxd', type=int, required=True,
                        help='set maximum depth to trim the whole tree into several smaller clades')
    parser.add_argument('-f', '--format', type=str, default='newick',
                        help='format of the phylogenetic tree (default: newick)')
    parser.add_argument('--mins', type=int, default=3,
                        help='minimum size (no. of nodes) of a trimmed clade. It must be >=3 (default: 3)')
    parser.add_argument('--font', type=int, default=6,
                        help='set font size for the size of text in the phylogenetic tree (default: 6)')
    parser.add_argument('-g', '--group', type=str, default='grouped_species',
                        help='output of each node in selected group (default: grouped_species)')
    parser.add_argument('--figure', type=str, default='tree.pdf',
                        help='name of figure for output where the selected clades are colored (default: tree.pdf)')
    args = parser.parse_args()
    biotree = Tree(args.tree, args.format)
    biotree.grouping(args.maxd, args.mins)
    print('Trimming the whole tree into smaller grpups...')
    sum_nodes = 0
    ofile = open(args.group, 'w')
    for idx, x in enumerate(biotree.grouped_nodes):
        value = '\t'.join(x)
        ofile.write('%s\n' % value)
        sum_nodes += len(x)
    print('Totally, %d species were used in %d groups' % (sum_nodes, idx+1))
    print('The grouped species has been written to %s' % args.group)
    biotree.save(args.figure, args.font)
    print('The tree has been shown in %s' % args.figure)

