# obtaining 1-to-1 orthologous relationships
#
# Chun-Ping Yu  chpngyu@gmail.com
# Dec 2, 2016
#


epilog_info = '''The 1-to-1 orthologous relationships of a group of OrthoFinder are determined by
there is only one gene in one species and
the species of these genes have to be included in all preprovided species (given by option --species).
'''

def assign(orthoships, all_species):
    _species2genes = {}
    for orth in orthoships:
        _species, _gene = orth.split('|')
        if _species in _species2genes:
            return False, None
        else:
            _species2genes[_species] = _gene
    common_species = all_species & set(_species2genes.keys())
    if len(common_species) == len(all_species):
        return True, _species2genes
    return False, None

import argparse
parser = argparse.ArgumentParser(description='Obtaining 1-to-1 orthologous relationships',
                                 epilog=epilog_info)
parser.add_argument('-o', '--orthology', type=str, required=True,
                    help='Orthologous groups from the output of OrthoFinder')
parser.add_argument('-s', '--species', type=str, required=True,
                    help='Provide a list of species be used. Only the gene group of orthologs included in all spepces will be reported (see below)')
parser.add_argument('-r', '--result', type=str, default='result.txt',
                    help='Result of output in a tab-delimited format with 3 columns (default: result.txt)')
args = parser.parse_args()

species = set([s.strip() for s in open(args.species) if s.strip() != ''])
print('There are %d speices in %s' % (len(species), args.species))

ofile = open(args.result, 'w')
ofile.write('Species\tGene\tGroup\n')
idx = 1
total = 0
with open(args.orthology, 'U') as ifile:
    for line in ifile:
        total += 1
        successful, species2genes = assign(line.split()[1:], species)
        if successful:
            for s in species2genes:
                ofile.write('%s\t%s\t%d\n' % (s, species2genes[s], idx))
            idx += 1
ofile.close()
print('%d(/%d in total) 1-to-1 orthologous groups were obtained' % (idx-1, total))
