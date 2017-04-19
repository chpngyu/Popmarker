# utility
#
# Chun-Ping Yu chongyu@gmail.com
# Nov. 22, 2016
#

import sys
from math import log
from scipy.stats import chi2

class Progress:
    def __init__(self, total_no_step):
        self.total_no_step = float(total_no_step)
        self.bar_width = 50
        self.current_step = 0
        self.previous_percentage = 0
        sys.stdout.write('0%   10   20   30   40   50   60   70   80   90   100%\n')
        sys.stdout.write('|----|----|----|----|----|----|----|----|----|----|\n')

    def restart(self, total_no_step):
        sys.stdout.flush()
        self.__init__(total_no_step)

    def walk(self, stride=1):
        self.current_step += stride
        percentage = int(self.bar_width * (self.current_step/self.total_no_step))
        if percentage > self.previous_percentage and percentage <= self.bar_width:
            for i in range(self.previous_percentage, percentage):
                sys.stdout.write('*')
            sys.stdout.flush()
            self.previous_percentage = percentage
            if percentage == self.bar_width:
                sys.stdout.write('*\n')
                sys.stdout.flush()


def fisher_combined_probability_test(p_values):
    if 0 in p_values:
        return 0
    df = 2*len(p_values)
    chi_squred = -2.0*sum([log(p) for p in p_values])
    return chi2.sf(chi_squred, df)

def tied_rank(vec):
    _rank = 1
    _ranks = [_rank]
    prev_v = vec[0]
    for v in vec[1:]:
        _rank += 1
        if v != prev_v:
            yield _ranks
            _ranks = [_rank,]
        else:
            _ranks.append(_rank)
        prev_v = v
    yield _ranks

def rank_data(vec):
    idx_data = sorted([(u, v) for u, v in enumerate(vec)], key=lambda x: x[1])
    idx = [u[0] for u in idx_data]
    sorted_data = [u[1] for u in idx_data]
    _ranks = []
    for rnk in tied_rank(sorted_data):
        rnk_len = len(rnk)
        if rnk_len > 1:
            tied = sum(rnk)/float(rnk_len)
            _ranks.extend([tied]*rnk_len)
        else:
            _ranks.append(rnk[0])
    ranks = [0]*len(vec)
    for i, _i in enumerate(idx):
        ranks[_i] = _ranks[i]
    return ranks


def read_group_block(file_path):
    ifile = open(file_path, 'U')
    # read header
    tokens = ifile.readline().strip().split('\t')
    if len(tokens) != 3:
        raise RuntimeError('%s format error.' % file_path)
    species, gene, group = ifile.readline().strip().split('\t')
    current_group = group
    gene_family = [(species, gene, group),]
    for line in ifile:
        line = line.strip()
        if len(line) == 0:
            continue
        species, gene, group = line.split('\t')
        if group != current_group:
            gene_family = sorted(gene_family)
            all_species = [s[0] for s in gene_family]
            genes = [s[1] for s in gene_family]
            gene_family = [(species, gene),]
            yield all_species, genes, current_group
            current_group = group
        else:
            gene_family.append((species, gene))
    ifile.close()
    gene_family = sorted(gene_family)
    all_species = [s[0] for s in gene_family]
    genes = [s[1] for s in gene_family]
    yield all_species, genes, group

from collections import OrderedDict
def read_orthologous_groups(file_path):
    group_block = read_group_block(file_path)
    species, genes, group = next(group_block)
    if len(species) == 0:
        raise RuntimeError('No data in %s' % file_path)
    idx = 0
    orthogene = OrderedDict()
    orthogene[group] = genes
    prev_group = group
    for next_species, genes, group in group_block:
        miss = set(next_species) - set(species)
        if len(miss):
            what_species = ', '.join(miss)
            msg = '{} species found in group {} but not in group {}: {}'.format(len(miss), group, prev_group, what_species)
            raise RuntimeError(msg)
        orthogene[group] = genes
        species = next_species
        prev_group = group
    return orthogene, species

from scipy.stats import pearsonr, spearmanr, kendalltau
def correlation_adaptor(x, y, method):
    '''select correlation functions from pcc, spearman, and kendall's tau.
        return correlation and single-tailed p-value (larger)'''
    if method == 0:
        corr, pvalue = pearsonr(x, y)
    elif method == 1:
        corr, pvalue = spearmanr(x, y)
    elif method == 2:
        corr, pvalue = kendalltau(x, y)
    else:
        raise ValueError('unrecognized method ' + str(method) + ', must be an integer from 1 to 4')
    if isnan(corr):
        corr = 0
        pvalue = 1
    elif corr < 0:
        pvalue = 1-pvalue*0.5
    else:
        pvalue *= 0.5
    return corr, pvalue

