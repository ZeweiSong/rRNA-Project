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
    while len(seq) > 0:
        frag = rnd_frag(minL, maxL) # generate a random length fragment
        pool.append(seq[:frag]) # add the fragment to pool
        seq = seq[frag:]
    return pool

# Generate the LFR sequencing data
# Assume that only a portion of sequence can be recovered (10-30%)
# I assume SE100
def lfr_data(pool, rate, readL):
    import random
    total_length = len(''.join(pool)) # the total length of the molecue
    print(total_length)
    expected_length = int(total_length * rate)
    expected_number = expected_length // readL
    print(expected_number)
    print(len(pool))
    all_read = [i[:readL] for i in pool] # All possible short reads
    return random.sample(all_read, expected_number)