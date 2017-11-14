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
from rrna_parser import rrna
from rrna_parser import data_io
import os

genome = 'GCF_000146045.2_R64'
fna_file = genome + '_genomic.fna'
gff_file = genome + '_genomic.gff'
fna_db_file = genome + '_genomic.fnadb'
gff_db_file = genome + '_genomic.gffdb'

if os.path.isfile(fna_db_file): # If the fnadb file is already there
    fna = data_io.read_fna_db(fna_db_file)
else: # If it is not there, creadte a new database, and write to a new fnadb file
    fna = data_io.create_fna_db(data_io.read_fasta_multiline(fna_file))
    data_io.write_fna_db(fna, fna_db_file)

if os.path.isfile(gff_db_file): # Check if gffdb exist
    gff = data_io.read_gff_db(gff_db_file)
else: # If not, read in the GFF file and create a new gff file
    gff = data_io.create_gff_db(gff_file, gff_db_file)