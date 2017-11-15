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


# Return the Xmer minimizer of the given sequence
def minimizer(seq, kmer=7):
    seq_rev = data_io.revcomp(seq)
    seqL = len(seq)
    min = 'T' * kmer # this is the last minimizer alphanetically

    for i in range(0, seqL - kmer + 1):
        current_kmer = seq[i:i+kmer] # Get the current kmer sequence
        if current_kmer < min:
            min = current_kmer
        current_kmer = seq_rev[i:i+kmer]
        if current_kmer < min:
            min = current_kmer
    return min


# A slider iterator that use a X size window to scan a sequence from 5' to 3'
class Slider:
    def __init__(self, seq, size):
        self.size = size
        self.seq = seq
        self.start = 0
        self.stop = self.start + size

    def __iter__(self):
        return self

    def __next__(self):
        self.start += 1
        self.stop += 1
        if self.stop > len(self.seq):
            raise StopIteration
        else:
            return self.seq[self.start:self.stop]