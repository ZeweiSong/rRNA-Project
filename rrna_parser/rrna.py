#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Fri Nov  3 16:29:34 2017

Zewei Song
BGI-Shenzhen
songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division
from rrna_parser import data_io
import json
import re
import random

class rrna_miner:
    def __init__(self, fnadb, gffdb): # Need fna_json, gffdb_gttutils, and the Class primer
        self.fnadb = fnadb
        self.gffdb = gffdb
    # Find all primers site on all contigs with max mismatch = 1
    def mine_primers(self, primer_class, mismatch = 1):
        loci = {}
        primer_list = list(primer_class)
        print(primer_list)
        for contig, seq in self.fnadb.items():
            print('Searching {0}: {1} bp ...'.format(contig, len(seq)))
            loci[contig] = []  # Create a new dict for each contig
            for p in primer_list:
                pos = data_io.search_oligo_loc(seq, p[1],
                                               mismatch=mismatch)  # Find all loci of current primer for current contig
                for key, value in pos.items():
                    for item in value:
                        if key.isupper():
                            strand = '+'
                        elif key.islower():
                            strand = '-'
                        else:
                            strand = 'm'
                        loci[contig].append([p[0], key, item[0], item[1], item[1] - item[0] + 1, strand, 'p'])

        count = []
        for key, value in loci.items():
            count.append([key, len(value)])
        count.sort(key=lambda x: x[1], reverse=True)
        print('Contig {0} has the most ({1}) primer site.'.format(count[0][0], count[0][1]))

        return loci

    # Find all rRNA genes on all contigs based on GFF annotation
    def mine_rrna_genes(self):
        rRNA = {}
        for gene in self.gffdb.features_of_type('gene'):  # Find all genes
            if gene.attributes['gene_biotype'][0] == 'rRNA':  # check if the gene is a rRNA
                for i in self.gffdb.children(gene, featuretype='rRNA'):
                    product = i.attributes['product']
                try:
                    rRNA[gene[0]].append(
                        [gene['Name'][0], product[0], gene.start, gene.end, gene.end - gene.start + 1, gene.strand,
                         'g'])
                except KeyError:
                    rRNA[gene[0]] = [
                        [gene['Name'][0], product[0], gene.start, gene.end, gene.end - gene.start + 1, gene.strand,
                         'g']]

        for key, value in rRNA.items():
            rRNA[key].sort(key=lambda x: x[2])
        return rRNA

    # Find the perfect rRNA with all expected primers and genes in one contig
    def mine(self, primers, rrna_genes):
        combined = {}
        for contig in list(self.fnadb.keys()):  # search for all contigs:
            try:
                combined[contig] = primers[contig] + rrna_genes[contig]
                combined[contig].sort(key=lambda x: x[2])
            except KeyError:
                combined[contig] = []

        # Identify rRNA using regular expression pattern
        # plus = 'gppgppgppp' or minus = 'gpppgppgpp'
        # Beware this is a very stringent method (may not work for poor annotated genomes)
        plus = 'gppgppgppp'
        minus = 'gpppgppgpp'
        rrna_candidate = {}

        for contig in self.fnadb.keys():
            rrna_plus = [i for i in combined[contig] if (i[-1] == 'g' and i[-2] == '+') or i[-1] == 'p']
            rrna_minus = [i for i in combined[contig] if (i[-1] == 'g' and i[-2] == '-') or i[-1] == 'p']

            for m in re.finditer(plus, ''.join([i[-1] for i in rrna_plus])):
                print(m)
                try:
                    rrna_candidate[contig].append(rrna_plus[m.start():m.end()])
                except KeyError:
                    rrna_candidate[contig] = [rrna_plus[m.start():m.end()]]
            for m in re.finditer(minus, ''.join([i[-1] for i in rrna_minus])):
                print(m)
                try:
                    rrna_candidate[contig].append(rrna_minus[m.start():m.end()])
                except KeyError:
                    rrna_candidate[contig] = [rrna_minus[m.start():m.end()]]
        print(rrna_candidate)

        # Get the locations all on the plus strand (as consistent with the input fna file)
        # Generate the rrna_predicted dictionary (can be save as a json file)
        rrna_predicted = {}  # it's element can be used as a rrna Class

        for contig, rrna_list in rrna_candidate.items():
            rrna_predicted[contig] = {}  # rrna list for current contig
            i = 1
            for rrna in rrna_list:  # iterate through all predicted rRNA section in current chromosome
                rrna_relocate = []
                rrna_name = contig + '_rrna' + str(i)
                i += 1

                # Get the start, stop, and strand information
                start = min([i[2] for i in rrna])  # find the start
                stop = max([i[3] for i in rrna])  # find the stop
                print(start, stop)
                strand = list(set([i[-2] for i in rrna if i[-1] == 'g']))[0]
                print(strand)

                # Change the relative position of all genes/primers in the same rRNA section
                # Also reverse compliment primers if the rrna is on minus strand
                for gene in rrna:
                    gene_relocate = []
                    if strand == '+':
                        for index, element in enumerate(gene):
                            if index == 2:  # This is the start of the original
                                gene_relocate.append(gene[2] - start + 1)
                            elif index == 3:  # This is the stop of the original
                                gene_relocate.append(gene[3] - start + 1)
                            else:
                                gene_relocate.append(element)

                    elif strand == '-':
                        start_new = abs(gene[3] - stop) + 1
                        stop_new = abs(gene[2] - stop) + 1
                        for index, element in enumerate(gene):
                            if index == 2:  # this is the start of the original gene/primer
                                gene_relocate.append(start_new)
                            elif index == 3:  # This is the stop of the original gene/primers
                                gene_relocate.append(stop_new)
                            elif index == 4:  # Calculate the new length, have to be the same as the old
                                gene_relocate.append(stop_new - start_new + 1)
                            elif index == 1 and gene[-1] == 'p':  # This is a primer
                                gene_relocate.append('-')
                            elif index == 5 and gene[-1] == 'p':  # This is a primer
                                if gene[5] == '+':
                                    gene_relocate.append('-')
                                elif gene[5] == '-':
                                    gene_relocate.append('+')
                            else:
                                gene_relocate.append(element)
                    rrna_relocate.append(gene_relocate)
                rrna_relocate.sort(key=lambda x: x[2])

                # Get the corresponding sequences on the original strand (or plus)
                print(start, stop, strand)
                if strand == '+':
                    seq = data_io.seq_slice(self.fnadb[contig], start, stop)
                elif strand == '-':
                    seq = data_io.seq_slice(self.fnadb[contig], start, stop, reverse_compliment=True)
                rrna_predicted[contig][rrna_name] = {'seq': seq, 'gff': rrna_relocate}
        return rrna_predicted

class primer:
    def __init__(self):
        self.primer_set = []
        self.index = 0

    def __iter__(self):
        return self

    def __next__(self):
        if len(self.primer_set) == 0:
            print('There is no primer, add some use add_primer().')
        else:
            try:
                current_primer = self.primer_set[self.index] # return a primer tuple
            except IndexError: # Stop at the last primer
                raise StopIteration
            self.index += 1
        return current_primer

    def add_primer(self, name, seq): # Add a new primer tuple
        self.primer_set.append((name, seq))
        
    def list(self):
        print('Primer_name\tSequence\n')
        for primer in self.primer_set:
            print('{0:s}\t{1:s}\n'.format(primer[0], primer[1]))
        return None

    def count(self):
        return(len(self.primer_set))


class ribosomeRNA:
    def __init__(self, rrna): # predicted_rrna is one rRNA section from the dictionary rrna_predicted[contig][rrna]
        self.seq = rrna['seq'] # the DNA sequence
        rrna_components = rrna['gff'] # This is a list of gene/primers, thus components
        self.gene = [i for i in rrna_components if i[-1] == 'g']
        self.primer = [i for i in rrna_components if i[-1] == 'p']
        self.primer_dict = {'f':{}, 'r':{}}
        for p in self.primer:
            if p[-2] == '+':
                self.primer_dict['f'][p[0]] = (p[2], p[3])
            elif p[-2] == '-':
                self.primer_dict['r'][p[0]] = (p[2], p[3])

    # return a list of all primers
    def primers(self):
        return [i[0] for i in self.primer]
    def forward_primers(self):
        return [i[0] for i in self.primer if i[-2] == '+']
    def reverse_primers(self):
        return [i[0] for i in self.primer if i[-2] == '-']
    def amplicon(self, f, r):
        from rrna_parser import data_io
        if f not in self.forward_primers():
            print('{0} is not in the froward primer pool.'.format(f))
        elif r not in self.reverse_primers():
            print('{0} is not in the reverse primer pool.'.format(r))
        else:
            start = self.primer_dict['f'][f][0]
            stop = self.primer_dict['r'][r][1]
        return data_io.seq_slice(self.seq, start, stop, reverse_compliment=False)

# from the predicted rrna list create a json file to save in .rrnadb
def create_rrna_db(rrna_predicted, json_file):
    with open(json_file, 'w') as f:
        json.dump(rrna_predicted, f)
    return None

# Return the minimizer of a given kmer sequence by length ml
def minimizer(kmer, ml):
    kmer
    kmer_rev = data_io.revcomp(kmer)

    mlist = []
    for i in range(0, len(kmer) - ml + 1):
        mlist.append(kmer[i:i + ml])
        mlist.append(kmer_rev[i:i + ml])

    mlist.sort()
    return mlist[0]


# Random mutate a given sequence with on the given number of nucleotide
# If with_replacement=True, the mutation happens one at a time, so multiulpe mutations can occur on one loci.
def mutation(seq, number, with_replacement=False):
    i = range(len(seq))
    seq = [i for i in seq]
    base = {'A':['T','C','G'],'T':['A','C','G'],\
            'C':['A','T','G'],'G':['A','T','C']}

    if with_replacement:
        loci = random.choices(i, number)
    else:
        loci = random.sample(i, number)

    for j in loci:
        seq[j] = random.sample(base[seq[j]], 1)[0]
    return ''.join(seq)