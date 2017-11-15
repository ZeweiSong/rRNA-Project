#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 17:15:34 2017

@author: Zewei Song (songzewei@outlook.com)
"""
#%%
from __future__ import print_function
from __future__ import division

# randomly generate
def rnd_frag(a, b):
    import random
    return random.randint(a, b)

# Shear the target DNA with random length in a range
def fragmentation(seq, minL, maxL):
    pool = []
    start = 1
    while len(seq) > 0:
        frag = rnd_frag(minL, maxL) # generate a random length fragment
        stop = start + frag - 1
        pool.append((seq[:frag], start, stop)) # add the fragment to pool
        start = stop + 1
        seq = seq[frag:]
    return pool

# Generate the LFR sequencing data
# Assume that only a portion of sequence can be recovered (10-30%)
# I assume SE100
def stlfr_data(pool, rate, readL):
    import random
    total_length = len(''.join([i[0] for i in pool])) # the total length of the molecue
    print(total_length)
    expected_length = int(total_length * rate)
    expected_number = expected_length // readL
    print(expected_number)
    print(len(pool))
    all_read = [(i[0][:readL], i[1], i[1] + readL - 1) for i in pool] # All possible short reads
    run = random.sample(all_read, expected_number)
    run.sort(key=lambda x:x[2])
    return run

# Write to a FASTA file
def write_stlfr(run, filename, species='Organism', magbead = 0):
    with open(filename, 'w') as f:
        for item in run:
            label = '>' + species + '-' + str(magbead) + '-' + str(item[1]) + '-' + str(item[2])
            f.write('{0}\n'.format(label))
            f.write('{0}\n'.format(item[0]))