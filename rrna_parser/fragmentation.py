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
# Setting minL > 100 will result in gap in the final data (will resolve later)
def fragmentation(seq, minL, maxL):
    pool = []
    start = 1
    while len(seq) >= minL:
        frag = rnd_frag(minL, maxL) # generate a random length fragment
        stop = start + frag - 1
        pool.append((seq[:frag], start, stop)) # add the fragment to pool
        start = stop + 1
        seq = seq[frag:]
    return pool

# Generate the LFR sequencing data
# Assume that only a portion of sequence can be recovered (10-30%)
# I assume SE100
def stlfr_data(pool, rate, readL, magbead = 0, sample = 'sample'):
    import random
    total_length = len(''.join([i[0] for i in pool])) # the total length of the molecue
    #print('The total length of template DNA is {0} bp.'.format(total_length))
    expected_length = int(total_length * rate)
    expected_number = expected_length // readL + 1
    #print('We want to recover {0}% = {1} bp sequences, using {2} bp read.'.format(rate*100.0, expected_length, readL))
    #print('Picking {0} reads from the total of {1} reads ...'.format(expected_number, len(pool)))
    all_read = [(i[0][:readL], i[1], i[1] + readL - 1, magbead) for i in pool] # All possible short reads, number of magbead is saved in 4th position
    run = random.sample(all_read, expected_number)
    run.sort(key=lambda x:x[2])
    
    # Convert the run data into the list/tuple format
    run_output = []
    for item in run:
        run_output.append((sample + '-' + str(magbead) + '-' + str(item[1]) + '-' + str(item[2]), item[0]))
    return tuple(run_output)


# Generate a sequence bin that simulate the data of stLFR (short read on a single molecular)
def sequence_bin(seq, minL, maxL, rate, readL):
    pass

# Write to a FASTA file
def write_stlfr(run, filename, species='Organism'):
    with open(filename, 'w') as f:
        for item in run:
            label = '>' + species + '-' + item[0]
            f.write('{0}\n'.format(label))
            f.write('{0}\n'.format(item[1]))