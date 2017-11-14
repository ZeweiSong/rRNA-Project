#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Tue Nov  7 17:10:41 2017

Zewei Song
BGI-Shenzhen
songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division

# Need to move all gff related content here:

# Try to find all rRNA in the GFF database
# The GGF may not formated as the one I tested
def rrna(db):
    rRNA_gene = {}
    sequence_region = {}
    for gene in db.features_of_type('gene'):
        if gene.attributes['gene_biotype'][0] == 'rRNA':
            rRNA_gene[gene['Name'][0]] = gene
            # Then we count the gene per contig
            try:
                sequence_region[gene[0]] += 1
            except KeyError:
                sequence_region[gene[0]] = 1

    for key, value in sequence_region.items():
        print('{0}: {1}'.format(key, value))
    return rRNA_gene
