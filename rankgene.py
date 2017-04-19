# ranking genes based on the correlation between phylogenetic tree and sequence similarities
#
# Chun-Ping Yu  chpngyu@gmail.com
# July 26, 2016
#
# revised and included p-value on Nov. 22, 2016
# included correlations in clades on Dec. 26, 2016
#

import numpy as np
np.seterr(divide='ignore', invalid='ignore')
import sys
from collections import defaultdict
from score import Orthoship
from biotree import Tree
from math import isnan
from random import shuffle
from utility import Progress, fisher_combined_probability_test, rank_data, correlation_adaptor

try:
    from statsmodels.sandbox.stats.multicomp import multipletests
    do_fdr = True
except ImportError:
    print('No StatsModels installed, the FDR estimation will be performed.')
    do_fdr = False

epilog_info = '''The available methods of correlation include
  0: Pearson correlation coefficient
  1: Spearman rank-order correlation coefficient
  2: Kendall's tau

The methods for adjustment of p-values (FDR) include
  0: bonferroni: one-step correction
  1: sidak: one-step correction
  2: holm-sidak: step down method using Sidak adjustments (default)
  3: holm: step-down method using Bonferroni adjustments
  4: simes-hochberg: step-up method  (independent)
  5: hommel: closed method based on Simes tests (non-negative)
  6: fdr_bh: Benjamini/Hochberg  (non-negative)
  7: fdr_by: Benjamini/Yekutieli (negative)
  8: fdr_tsbh: two stage fdr correction (non-negative)
  9: fdr_tsbky: two stage fdr correction (non-negative)'''

import argparse
parser = argparse.ArgumentParser(description="Ranking genes by correlations between sequences similarities and a reference tree",
                                 epilog=epilog_info, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-b', '--blast', type=str, required=True,
                    help='results of blastall in tabular format (-outfmt 6)')
parser.add_argument('-t', '--tree', type=str, required=True,
                    help='reference tree')
parser.add_argument('-f', '--format', type=str, default='newick',
                    help='format of the reference tree (default: newick)')
parser.add_argument('-o', '--ortholog', type=str, required=True,
                    help='relationships of orthology among species')
parser.add_argument('--clade', type=str,
                    help='also calculate correlations of species in the given clade(s)')
parser.add_argument('-r', '--result', type=str, default='result.txt',
                    help='output of result (default: result.txt)')
parser.add_argument('--corr', type=int, default=0, choices=range(3),
                    help='method for computing correlation. The default is pearson correlation(0); see below')
parser.add_argument('--fdrm', type=int, default=6, choices=range(10),
                    help='method used for adjustment of p-values (default is 6; see below)')
args = parser.parse_args()

biotree = Tree(args.tree, args.format)
if biotree.tree.count_terminals() < 3:
    raise RuntimeError('The no. of node in the tree is too small (<3)')

ortho = Orthoship(args.blast, args.ortholog)
if len(ortho.ortholog) == 0:
    raise RuntimeError('No any species found in the %s' % args.blast) 
first_species = next(iter(ortho.ortholog.keys()))

genes = ortho.ortholog[first_species].keys()
all_species = set(ortho.species.keys()).intersection([x.name for x in biotree.tree.get_terminals()])
print('There are %d orthologous gene pairs among %d species' % (len(genes), len(all_species)))

if args.clade:
    tmp_species_in_clade = defaultdict(set)
    species_in_clade = {}
    with open(args.clade, 'U') as ifile:
        for line in ifile:
            line = line.strip()
            if line == '':
                continue
            tokens = line.split('\t')
            tmp_species_in_clade[tokens[1]].add(tokens[0])
    # check whether species are in the leaves of the given tree
    _all_speices_in_clade = set()
    _no_clade = len(tmp_species_in_clade)
    _no_used_species = 0
    for clade in tmp_species_in_clade:
        _all_speices_in_clade.update(tmp_species_in_clade[clade])
        _used_species = tmp_species_in_clade[clade] & all_species
        if len(_used_species) > 2:
            species_in_clade[clade] = list(_used_species)
            _no_used_species += len(_used_species)
        else:
            print('No or too few species (must be >2) in %s clade. Got only %d species' % (clade, len(_used_species)))
    del tmp_species_in_clade
    if len(species_in_clade):
        print('From %s, there were %d species in %d clade(s), where %d species in %d clade(s) were used'
              % (args.clade, len(_all_speices_in_clade), _no_clade, _no_used_species, len(species_in_clade)))
        do_clade = True
else:
    print('The calculation of correlations within clade(s) will not be performed')
    do_clade = False

all_species = list(all_species) # tarnsfrom to list for later use
tree_distance = biotree.distance(all_species)

# computing the correlations between the reference tree and sequence similarities of genes
result = []
corr_method = ['Pearson correlation', 'Spearman correlation', "Kendall's tau"][args.corr]
print('Computing correlations (%s) and p-values for whole tree' % corr_method)
prgs = Progress(len(genes))
for gene in genes:
    score = ortho.get_scores_from(gene, all_species)
    pcc, pvalue = correlation_adaptor(tree_distance, score, args.corr)
    orthologous_genes = []
    for species in all_species:
        orthologous_genes.append(ortho.ortholog[species][gene])
    result.append((orthologous_genes, pcc, pvalue))
    prgs.walk()

if do_clade:
    print('Computing correlations (%s) and p-values for %d clade(s)...' % (corr_method, len(species_in_clade)))
    clade_tree_distance = []
    for clade in species_in_clade:
        _species = species_in_clade[clade]
        clade_tree_distance.append(biotree.distance(_species))
    prgs.restart(len(genes)*len(species_in_clade))
    for gene_idx, gene in enumerate(genes):
        pcc_by_clade = []
        pvalue_by_clade = []
        for clade_idx, clade in enumerate(species_in_clade):
            _species = species_in_clade[clade]
            score = ortho.get_scores_from(gene, _species)
            pcc, pvalue = correlation_adaptor(clade_tree_distance[clade_idx], score, args.corr)
            pcc_by_clade.append(pcc)
            pvalue_by_clade.append(pvalue)
            prgs.walk()
        avg_pcc = sum(pcc_by_clade)/len(pcc_by_clade)
        cmb_pvalue = fisher_combined_probability_test(pvalue_by_clade)
        result[gene_idx] += (avg_pcc, cmb_pvalue)
    clade_rank = rank_data([(r[4], -r[3]) for r in result])

pvalues = [r[2] for r in result]
# sorting by p-value and -pcc
rank = rank_data([(r[2], -r[1]) for r in result])
if do_fdr:
    rank_idx = len(all_species)+3
else:
    rank_idx = len(all_species)+2

if do_fdr:
    sys.stdout.write('Computing FDRs....')
    # computing FDRs
    fdr_method = ['bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel',
                  'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky'][args.fdrm]
    fdr = multipletests(pvalues, method=fdr_method)[1]
    if do_clade:
        pvalues = [r[4] for r in result]
        clade_fdr = multipletests(pvalues, method=fdr_method)[1]
        clade_rank = rank_data([(r[4], -r[3]) for r in result])
    sys.stdout.write('done\n')

# combining all results into a list, and
# sorting genes by the ranks (or average ranks if the clade calculation is set)
tmp_result = []
if do_clade:
    if do_fdr:
        for a, b, c, d, e in zip(result, fdr, rank, clade_fdr, clade_rank):
            tmp_result.append( (tuple(a[0]) + a[1:3] + (b, c) + a[-2:] + (d, e)) )
    else:
        for a, b, c in zip(result, rank, clade_rank):
            tmp_result.append( (tuple(a[0]) + a[1:3] + (b,) + a[-2:] + (c,)) )
else:
    if do_fdr:
        for a, b, c in zip(result, fdr, rank):
            tmp_result.append( (tuple(a[0]) + a[1:3] + (b, c)) )
    else:
        for a, b in zip(result, rank):
            tmp_result.append( (tuple(a[0]) + a[1:3] + (b,)) )
result = sorted(tmp_result, key=lambda x: x[rank_idx]) # sorting by rank
del tmp_result

# writing results to the putput
outfile = open(args.result, 'w')
head = '\t'.join(all_species)
if do_fdr:
    outfile.write('ID\t%s\t%s\tP-value\tFDR\tRank' % (head, corr_method))
    if do_clade:
        outfile.write('\t%s in clade\tCombined p-value in clade\tFDR in clade\tRank in clade' % corr_method)
else:
    outfile.write('ID\t%s\t%s\tP-value\tRank' % (head, corr_method))
    if do_clade:
        outfile.write('\t%s in clade\tCombined p-value in clade\tRank in clade' % corr_method)
outfile.write('\n')

for idx, r in enumerate(result):
    value_str = '\t'.join([str(u) for u in r])
    rank_id = 'ID{:>06d}'.format(idx+1)
    outfile.write('%s\t%s\n' % (rank_id, value_str))
outfile.close()
print('The results have been written to %s' % args.result)
