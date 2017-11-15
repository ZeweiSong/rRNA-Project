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
import importlib # This is used to reload module during testing, see Evernote for more details

# Read in the FASTA and GFF data.
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
    gff = data_io.read_gff_db(gff_db_file, report=True)
else: # If not, read in the GFF file and create a new gff file
    gff = data_io.create_gff_db(gff_file, gff_db_file)

# Read in primers
primers = rrna.primer() # create a primer object
primer_fna = data_io.read_fasta('rRNA_primers.fasta') # Read in the FASTA file of primers
for record in primer_fna: # Add primers to the primer object
    primers.add_primer(record[0], record[1])

# Fina all locations of primers in the genome
#importlib.reload(rrna)
rrna_mine = rrna.rrna_miner(fna, gff) # Create the miner object
primer_loci = rrna_mine.mine_primers(primers, mismatch=1) # Find all possible location for all primers on all contigs
gene_loci = rrna_mine.mine_rrna_genes() # Find all possible genes annotated as rRNA (may NOT exist in many genomeS)
rrna_candidate = rrna_mine.mine(primer_loci, gene_loci)
#%% Save the rrna_candidate to json
data_io.write_json(rrna_candidate, genome + '_rrna.rrnadb')

#%% Next time start from here
from rrna_parser import rrna
from rrna_parser import data_io
from rrna_parser import fragmentation
import os
import importlib # This is used to reload module during testing, see Evernote for more details

genome = 'GCF_000146045.2_R64'
rrna_candidate_db = genome + '_rrna.rrnadb'
rrna_candidate = data_io.read_json(rrna_candidate_db)
for contig, rrna_list in rrna_candidate.items():
    print(contig)
    for r in rrna_list.keys():
        print(r)

sequence = rrna.ribosomeRNA(rrna_candidate['NC_001144.5']['NC_001144.5_rrna1'])
amplicon = sequence.amplicon(sequence.primers()[0], sequence.primers()[6]).upper()

#%% Generate stLFR simulated data
importlib.reload(fragmentation)
shear_dna = fragmentation.fragmentation(amplicon, 120, 200)
run = fragmentation.stlfr_data(shear_dna, 0.3, 100)
fragmentation
fragmentation.write_stlfr(run, genome + '_run1.fa', species='SC', magbead=0)