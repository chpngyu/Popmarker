import sys
no_record = 50000
class Orthoship:
    def __init__(self, blast_file, ortho_group_file, orth_group_type=None):
        if orth_group_type:
            self.read_orthologous_group(ortho_group_file, orth_group_type)
        else:
            self.read_orthologous_pairs(ortho_group_file)
        self.read_blast_scores(blast_file)
        self.normalize()

    def read_orthologous_pairs(self, file_path):
        # self.species is {'speciesA': a set of genes, 'speciesB': a set of genes, 'speciesC': a set of genes,...}
        self.species = {}
        # self.ortholog is {'speciesA': {homogroup 1: gene, ...}, 'speciesB': {homogroup 1: gene, ...}, ...} 
        self.ortholog = {}
        with open(file_path, 'U') as infile:
            infile.readline()
            for line in infile:
                # tokens[0], tokens[1], tokens[2] = species, gene, group
                tokens = line.strip().split('\t')
                if len(tokens) == 0:
                    continue
                elif len(tokens) != 3:
                    raise RuntimeError('format error in %s, it has to be species \tab gene \tab group' % file_path)
                if tokens[0] not in self.species:
                    self.species[tokens[0]] = set()
                    self.ortholog[tokens[0]] = {}
                self.species[tokens[0]].add(tokens[1])
                self.ortholog[tokens[0]][int(tokens[2])] = tokens[1]

    def read_orthologous_groups(self, file_path, orth_group_type):
        # self.species is {'speciesA': a set of genes, 'speciesB': a set of genes, 'speciesC': a set of genes,...}
        self.species = {}
        # self.ortholog is {'speciesA': {homogroup 1: gene, ...}, 'speciesB': {homogroup 1: gene, ...}, ...} 
        self.ortholog = {}
        with open(file_path, 'U') as infile:
            columns = infile.readline().strip().split()
            idx = columns.index('P-value')-1
            for sp in columns[1:idx]:
                self.species[sp] = set()
                self.ortholog[sp] = {}
            no_count = 0
            for line in infile:
                tokens = line.strip().split('\t')
                if len(tokens) == 0:
                    continue
                group = int(tokens[0][2:])
                for i, gene in enumerate(tokens[1:idx]):
                    sp = column[i+1]
                    self.species[sp].add(gene)
                    self.ortholog[sp][group] = gene
                no_count += 1
                if no_count > orth_group_type:
                    break
        

    def is_ortholgous(self, first_gene, second_gene):
        species, gene = first_gene.split('|')
        if (species not in self.species) or (gene not in self.species[species]):
            return False
        species, gene = second_gene.split('|')
        if (species not in self.species) or (gene not in self.species[species]):
            return False
        return True

    def read_blast_scores(self, file_path):
        self.scores = {}
        total_records = 0
        taken_records = 0
        with open(file_path, 'U') as infile:
            sys.stdout.write('Opening %s (*=%d records): ' % (file_path, no_record))
            for line in infile:
                total_records += 1
                tokens = line.strip().split('\t')
                if len(tokens) != 12: # not blast tabular format
                    continue
                if self.is_ortholgous(tokens[0], tokens[1]):
                    taken_records += 1
                    if (taken_records % no_record) == 0:
                        sys.stdout.write('*')
                    key = tuple(sorted([tokens[0], tokens[1]]))
                    evalue, bit_score = float(tokens[10]), float(tokens[11])
                    self.scores[key] = [evalue, bit_score]
            sys.stdout.write(' done\n')
        print('There are %d records in %s, where %d records are used' % (total_records, file_path, taken_records))

    def normalize(self):
        print('Scaling blast scores...')
        self.set_evalue_score()
        first_key = next(iter(self.scores))
        min_b = self.scores[first_key][1]
        max_b = min_b
        for key in self.scores:
            if self.scores[key][1] < min_b:
                min_b = self.scores[key][1]
            elif self.scores[key][1] > max_b:
                max_b = self.scores[key][1]
        for key in self.scores:
            s1 = self.scores[key][0]
            s2 = (self.scores[key][1]-max_b)/(min_b-max_b)
            self.scores[key] = 0.5*(s1+s2)
        print('done\n')

    def set_evalue_score(self):
        from math import log
        evalues = [self.scores[x][0] for x in self.scores if self.scores[x][0] != 0]
        if len(evalues) == 0:
            raise ValueError('All blast E-values are 0')
        min_e = evalues[0]
        max_e = evalues[0]
        for e in evalues:
            if e < min_e:
                min_e = e
            elif e > max_e:
                max_e = e
        print('Minimum non-zero E-value is %g\nMaximum E-value is %g' % (min_e, max_e))
        max_evalue = -log(min_e*0.1)
        min_evalue = -log(max_e)
        for key in self.scores:
            evalue = self.scores[key][0]
            if evalue != 0:
                self.scores[key][0] = (-log(evalue)-max_evalue)/(min_evalue - max_evalue)
            else:
                self.scores[key][0] = 0

    def get_scores_from(self, a_homologous_group, among_species):
        _scores = []
        for idx, speciesA in enumerate(among_species):
            geneA = speciesA + '|' + self.ortholog[speciesA][a_homologous_group]
            for speciesB in among_species[idx+1:]:
                geneB = speciesB + '|' + self.ortholog[speciesB][a_homologous_group]
                key = tuple(sorted([geneA, geneB]))
                if key in self.scores:
                    _scores.append(self.scores[key])
                else:
                    _scores.append(0.9)
        return _scores
