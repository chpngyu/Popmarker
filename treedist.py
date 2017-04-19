# computing distances (1-PCC) between trees
#
# Chun-Ping Yu  chpngyu@gmail.com
# Dec 3, 2016
#
from scipy.stats import pearsonr
from biotree import Tree
import argparse
parser = argparse.ArgumentParser(description='computing distances between trees')
parser.add_argument('-r', '--reference', type=str, required=True,
                    help='phylogenetic tree as reference')
parser.add_argument('-q', '--query', type=str, nargs='+',
                    help='other tree(s) for comparison(s)')
parser.add_argument('-f', '--format', type=str, default='newick',
                    help='format of the phylogenetic tree (default: newick)')
parser.add_argument('-t', '--targets', type=str, nargs='*', default=None,
                    help='only compute the distance for the giveb targets (must > 2 targets)')
parser.add_argument('-o', '--output', type=str, default='tree_distance.txt',
                    help='result ouptput for distances (default: tree_distance.txt)')
args = parser.parse_args()

ref_tree = Tree(args.reference, args.format)
ref_species = [x.name for x in ref_tree.tree.get_terminals()]
if len(ref_species) < 3:
    raise RuntimeError('The no. of node in the tree is too small (<3)')

if args.targets != None:
    if len(args.targets) > 2:
        ref_species = set(args.targets) & set(ref_species)
        check = set(args.targets) - set(ref_species)
        if len(check):
            missing_species = ', '.join(check)
            print('%d nodes are missing in %s: %s' %
                  (len(check), args.reference, missing_species))
        ref_species = list(ref_species)
    else:
        raise RuntimeError('--targets must give more than 2 targets')
ref_distance = ref_tree.distance(ref_species)

ofile = open(args.output, 'w')
ofile.write('Query\tDistance(1-PCC)\n')
for query in args.query:
    query_tree = Tree(query, args.format)
    if args.targets != None:
        query_species = set(args.targets) & set([x.name for x in query_tree.tree.get_terminals()])
    else:
        query_species = set([x.name for x in query_tree.tree.get_terminals()])
    print('Using %d target nodes to compute pairwise distances' % len(query_species))
    check = set(ref_species) - query_species
    if len(check):
        missing_species = ', '.join(check)
        print('%d nodes are missing in %s: %s' % (len(check), query, missing_species))
    else:
        query_distance = query_tree.distance(ref_species)
        pcc, pvalue = pearsonr(ref_distance, query_distance)
        ofile.write('%s\t%g\n' % (query, 1-pcc))
ofile.close()
