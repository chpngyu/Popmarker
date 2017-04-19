# computing correlations between combined genes and a reference tree
#
# Chun-Ping Yu  chpngyu@gmail.com
# Jan 7, 2017
#

from itertools import combinations
import argparse
parser = argparse.ArgumentParser(description="Computing correlations between combined genes and a reference tree")
parser.add_argument('-b', '--blast', type=str, required=True,
                    help='results of blastall in tabular format (-outfmt 6)')
parser.add_argument('-t', '--tree', type=str, required=True,
                    help='reference tree')
parser.add_argument('-f', '--format', type=str, default='newick',
                    help='format of the reference tree (default: newick)')
parser.add_argument('--maxn', type=int, default=3,
                    help='maximum number of genes to be combined (default: 3)')
parser.add_argument('--clade', type=str,
                    help='also calculate correlations of species in the given clade(s)')
parser.add_argument('-r', '--result', type=str, default='result.txt',
                    help='output of result (default: result.txt)')
parser.add_argument('--corr', type=int, default=0, choices=range(3),
                    help='method for computing correlation. The default is pearson correlation(0); see below')
parser.add_argument('--min', type=float, default=0.5,
                    help='only combined genes which have the correlation larger than this value will be reported (default: 0.5)')
args = parser.parse_args()

biotree = Tree(args.tree, args.format)
if biotree.tree.count_terminals() < 3:
    raise RuntimeError('The no. of node in the tree is too small (<3)')

ortho = Orthoship(args.blast, args.ortholog, 1)
if len(ortho.ortholog) == 0:
    raise RuntimeError('No any species found in the %s' % args.blast) 
first_species = next(iter(ortho.ortholog.keys()))

genes = ortho.ortholog[first_species].keys()
all_species = set(ortho.species.keys()).intersection([x.name for x in biotree.tree.get_terminals()])
print('There are %d orthologous gene pairs among %d species' % (len(genes), len(all_species)))
