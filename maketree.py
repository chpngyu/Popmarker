# constructing gene tree by giving gene sequences
#
# Chun-Ping Yu  chpngyu@gmail.com
# December 15, 2016
#

from subprocess import call
def run_mafft_commandline(_args, in_file, out_file):
    _cmd = ['mafft',]
    _cmd += ['--op', str(_args.op)]
    _cmd += ['--ep', str(_args.ep)]
    _cmd += ['--maxiterate', str(_args.maxiterate)]
    if _args.reorder:
        if _args.con or args.orthogrp:
            print('Since --con or --orthogrp has set, --reorder will not be used')
        else:
            _cmd.append('--reorder')
    if not _args.report:
        _cmd.append('--quiet')
    if _args.thread != None:
        _cmd += ['--thread', str(_args.thread)]
    if _args.accuracy != None:
        _cmd.append('--'+_args.accuracy)
    _cmd += [in_file, '>', out_file]
    _cmd = ' '.join(_cmd)
    err = call(_cmd, shell=True)
    if err:
        print('Failed to run mafft by commandline: %s' % _cmd)
    return err

def run_fasttree_commandline(_args, in_file, out_file):
    _cmd = ['FastTree',]
    if _args.wag:
        _cmd.append('-wag')
    if _args.gamma:
        _cmd.append('-gamma')
    if not _args.report:
        _cmd.append('-quiet')
    if _args.nt:
        if _args.gtr:
            _cmd += ['-gtr', '-nt', in_file, '>', out_file]
        else:
            _cmd += ['-nt', in_file, '>', out_file]
    else:
        _cmd += [in_file, '>', out_file]
    _cmd = ' '.join(_cmd)
    err = call(_cmd, shell=True)
    if err:
        print('Failed to run FastTree by commandline: %s' % _cmd)
    return err

from Bio import SeqIO # for sequence
def gaps_in_alignment(seq_file):
    total_length = 0
    no_gap = 0
    for record in SeqIO.parse(seq_file, 'fasta'):
        total_length += len(record)
        no_gap += record.seq.count('-')
    return 100.0*no_gap/total_length

def replace_genename_by_species():
    pass

from collections import OrderedDict
from os import remove
from utility import read_orthologous_groups
import argparse
parser = argparse.ArgumentParser(description='Constructing a tree by a gene sequence or more')
parser.add_argument('-g', '--grank', type=str,
                    help='a table of ranking genes associated with a phylogenetic tree (the output of rankgene.py)')
parser.add_argument('--orthogrp', type=str,
                    help='or give all orthologous genes from the output of orthoships.py. --con will be set')
parser.add_argument('-s', '--sequence', type=str, required=True,
                    help="a mapping table for species and it sequence's filename")
parser.add_argument('-t', '--top', type=int,
                    help='perform a gene tree for the concatenated top N sequences of genes')
parser.add_argument('-n', '--nth', type=int, nargs='+',
                    help='perform gene tree(s) for the given top n-th gene(s)')
parser.add_argument('--gap', type=float, default=20,
                    help='maximum percentage of gap extension to be included in the alignments (default: 20)')
parser.add_argument('--con', action='store_true',
                    help='concatenate the selected sequences (after alignment) into one sequence and make one tree')

parser.add_argument('--maxiterate', type=int, default=1000,
                    help='[Mafft] maximum number of iterative refinement (default: 1000)')
parser.add_argument('--op', type=float, default=1.53,
                    help='[Mafft] opening gap penalty (default: 1.53)')
parser.add_argument('--ep', type=float, default=0.0,
                    help='[Mafft] offset (works like gap extension penalty) (default: 0.0)')
parser.add_argument('--reorder', action='store_true',
                    help='[Mafft] output order (default: not set, use input order)')
parser.add_argument('--accuracy', type=str, nargs='?', const=None,
                    choices=['localpair', 'genafpair', 'globalpair'],
                    help='[Mafft] perform more accurate aligmmnent')
parser.add_argument('--thread', type=int, nargs='?', const=None,
                    help='[Mafft] number of threads (default: 1)')
parser.add_argument('--wag', action='store_true',
                    help='[FastTree] Whelan-And-Goldman 2001 model (amino acid alignments only)')
parser.add_argument('--gamma', action='store_true',
                    help='[FastTree] after optimizing the tree under the CAT approximation, rescale the lengths to optimize the Gamma20 likelihood')
parser.add_argument('--nt', action='store_true',
                    help='[FastTree] input sequence is nucleotide alignment')
parser.add_argument('--gtr', action='store_true',
                    help='[FastTree] generalized time-reversible model (nucleotide alignments only)')
parser.add_argument('--report', action='store_true',
                    help='[Mafft & FastTree] show detailed information of prgress (default: not show)')

args = parser.parse_args()

#### A. obtaining orthologous genes from the output of rankgene.py or orthoships.py
orthogene = OrderedDict()
if args.grank:
    if (args.top == None) and (args.nth == None):
        print('give option(s) for --top and/or --nth')
        print('Program terminated')
        exit(1)
    used_records = []
    if args.top:
        used_records += list(range(1, args.top+1))
    if args.nth:
        used_records += args.nth
    max_record = max(used_records)
    ifile = open(args.grank, 'U')
    tokens = ifile.readline().strip().split('\t')
    try:
        idx = tokens.index('P-value')
        last_item = idx-1
        species = tokens[1:last_item]
    except ValueError:
        print('Format error in %s (no item for P-value in the header)' % args.grank)
        print('Program terminated')
        exit(1)
    idx = 0
    for line in ifile:
        line = line.strip()
        if line == '':
            continue
        tokens = line.split('\t')
        idx += 1
        if idx in used_records:
            orthogene[tokens[0]] = tokens[1:last_item]
        elif idx > max_record:
            break
    ifile.close()
    print('Got %d group(s) of orthologous genes within %d species from %s' % (len(orthogene), last_item-1, args.grank))
elif args.orthogrp:
    orthogene, species = read_orthologous_groups(args.orthogrp)
    print('Got %d group(s) of orthologous genes within %d species from %s' % (len(orthogene), len(species), args.orthogrp))
else:
    print('give an option for either --grank or --orthogrp')
    print('Program terminated')
    exit(1)

#### B. obtaining sequences for the orthologous genes
seq_files = dict([tuple(line.strip().split('\t')) for line in open(args.sequence) if len(line.strip()) != 0])
print("Got %d pairs of species and sequenence filenames" % len(seq_files))

missing_species = set(species) - set(seq_files.keys())
if len(missing_species):
    print('Following sequences of species (%d) are missing:' % len(missing_species))
    print('%s' % '\n'.join(missing_species))
    print('Program terminated!')
    exit(1)

orthoseq = {}
for idx, sp in enumerate(species):
    genes = set([sp + '|' + orthogene[group][idx] for group in orthogene])
    seq_records = [(record.id, record) for record in SeqIO.parse(seq_files[sp], 'fasta') if record.id in genes]
    missing_gene = genes - set([s[0] for s in seq_records])
    if len(missing_gene):
        print('The following genes are missing (%d):' % len(missing_gene))
        print('%s' % '\n'.join(missing_species))
        print('Program terminated!')
        exit(1)
    print('Got %d sequences from %s' % (len(seq_records), seq_files[sp]))
    orthoseq[sp] = dict(seq_records)
        
#### C. performing multiple sequences alignment for each group of orthologous sequences
used_group = []
for idx, group in enumerate(orthogene):
    idx += 1
    print('Processing %s (%d/%d)' % (group, idx, len(orthogene)))
    seq_records = []
    for idx, gene in enumerate(orthogene[group]):
        sp = species[idx]
        gene_name = sp + '|' + gene
        record = orthoseq[species[idx]][gene_name]
        # change sequence name to species name
        record.id = sp
        seq_records.append(record)
    seq_outf = group + '.fasta'
    SeqIO.write(seq_records, seq_outf, 'fasta')
    # perform mafft
    aln_out = group + '.aln'
    run_mafft_commandline(args, seq_outf, aln_out)
    remove(seq_outf)
    # check gap percentage
    gap_perc = gaps_in_alignment(aln_out)
    if gap_perc < args.gap:
        used_group.append(aln_out)
    else:
        print('{0} was discarded due to higher gaps {1:5.2f}% (>= {2:5.2f}%)'.format(aln_out, gap_perc, args.gap))
        remove(aln_out)

if (args.con and (len(used_group) > 1)) or args.orthogrp:
    # concatenate all single gene alignments
    conseq_file = 'concatenated_seqs_%d.aln' % len(used_group)
    print('Concatenating %d alignments into one file (%s)' % (len(used_group), conseq_file))
    seqs = SeqIO.to_dict(SeqIO.parse(used_group[0], 'fasta'))
    remove(used_group[0])
    for seq_file in used_group[1:]:
        other_seqs = SeqIO.to_dict(SeqIO.parse(seq_file, 'fasta'))
        for seq_id in seqs:
            seqs[seq_id].seq += other_seqs[seq_id].seq
        remove(seq_file)
    SeqIO.write(seqs.values(), conseq_file, 'fasta')
    used_group = [conseq_file,]

#### D. using the alignment(s) to make tree(s)
for idx, group in enumerate(used_group):
    out_file = group.split('.')[0] + '.tree'
    idx += 1
    print('Making a tree to %s (%d/%d)' % (out_file, idx, len(used_group)))
    run_fasttree_commandline(args, group, out_file)
