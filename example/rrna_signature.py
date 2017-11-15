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
from rrna_parser import rrna
from rrna_parser import kmer
from rrna_parser import fragmentation
import importlib

genome = 'GCF_000146045.2_R64'
rrna_candidate_db = genome + '_rrna.rrnadb'
rrna_candidate = data_io.read_json(rrna_candidate_db)
for contig, rrna_list in rrna_candidate.items():
    print(contig)
    for r in rrna_list.keys():
        print(r)

sequence = rrna.ribosomeRNA(rrna_candidate['NC_001144.5']['NC_001144.5_rrna1'])
amplicon = sequence.amplicon(sequence.primers()[0], sequence.primers()[6]).upper()

#%% Scan the full length Sc rRNA for minimizer
summary = [['kmer_size', 'Streak', 'Unique']]
for size in [5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35]:
    seq_scan = kmer.Slider(sequence.seq.upper(), 100)
    min_list = []
    kmer_size = size
    report_file = 'sc_rrna_minimizer_kmer' + str(kmer_size) + '.txt'
    for i in seq_scan:
        min_list.append(kmer.minimizer(i, kmer=kmer_size))

    # Mark the streak for minimizer
    # Three columns: sequence, same as one above (0 yes, 1 no), same as another streak (0 yes, 1 no).
    min_report = [(min_list[0], 1, 1)]
    min_unique = [min_list[0]]
    a = 0
    b = 0
    for min in min_list[1:]:
        if min == min_unique[-1]: # same as last minimizer, continue of current streak
            a = 0
            b = 0
        else: # Start the count of a new streak
            a = 1
            if min in min_unique: # same as a previous unique minimizer (i.e. remote redundancy)
                b = 0 # current steak as occured previously
            else:
                b = 1
            min_unique.append(min)
        min_report.append((min, a, b))

    with open(report_file, 'w') as f:
        f.write('Minimizer\tSame as above\tRemote redundant\n')
        for line in min_report:
            f.write('{0}\t{1}\t{2}\n'.format(line[0], line[1], line[2]))
    streak = sum([i[1] for i in min_report])
    unique = sum([i[2] for i in min_report])
    summary.append((size, streak, unique))

with open('summary.txt', 'w') as f:
    for line in summary:
        f.write('%s\n' % '\t'.join([str(i) for i in line]))

#%% Generate stLFR simulated data
importlib.reload(fragmentation)
template = sequence.seq.upper()
for i in range(10):
    output = genome + '_run' + str(i) + '.fa'
    shear_dna = fragmentation.fragmentation(template, 120, 200)
    run = fragmentation.stlfr_data(shear_dna, 0.3, 100)
    fragmentation.write_stlfr(run, output, species='SC', magbead=i)
