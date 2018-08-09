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

# Read the entire content of the file by lines and store in list
# It is okay for reading a moderate size genome sequence file or a single profiling sample,
# but may not suitable for extremly large sequence file if it exceed the size of your memory.
# I have planned to write another module to read in file by trunk,
# but it is more tricky if FASTA file comes in multiple lines (which I believed is a legacy from Sanger age, and should be abandon)
def read_file(input_filename):
    with open(input_filename, 'rU') as f:
        content = f.readlines()  # f.read().split('\n') take more memory than f.readlines()
    return content # Return a list with all lines in the file. Lines contain \n at the end

# Read in a sequence file (FASTA or FASTQ) and store records in a list
def read_seqs(input_filename, file_type='default', output='default'):
    import sys
    # Check for record header
    with open(input_filename, 'rU') as f:
        content = f.readlines()
        head_symbol = content[0][0]

    if file_type == 'default':
        # Try to guess file format if not provided
        if head_symbol == '>':
            file_type = 'fasta'
        elif head_symbol == '@':
            file_type = 'fastq'
        else:
            print('%s is not a correct header for FASTA or FASTQ, please check you file.' % head_symbol)
            sys.exit()
    if file_type in ['fasta', 'Fasta', 'FASTA']:
        if content[0][0] == '@':
            print('%s seems to be a FASTQ file. Please use the correct format.' % input_filename)
            sys.exit()
        line_num = 2
        seq_line_num = 2
    elif file_type in ['fastq', 'Fastq', 'FASTQ']:
        if content[0][0] == '>':
            print('%s seems to be a FASTA file. Please use the correct format.' % input_filename)
            sys.exit()
        line_num = 4
        seq_line_num = 4
    else:
        print('Please specify a right format [fasta | fastq].')
        sys.exit()

    # set line number to 2 if want to export FASTA for a FASTQ file.
    if output in ['fasta', 'Fasta', 'FASTA']:
        seq_line_num = 2

    # Store record in list
    seq_num = int(len(content) / line_num)  # Total number of records
    output_content = []
    for i in range(seq_num):
        temp = []  # Store current record
        for line in range(seq_line_num):
            # Loop through each line in the current record (2 for fasta, 4 for fastq)
            temp.append(content[i * line_num + line][:-1])  # [:-1] remove '\n' at the end of each line
            content[i * line_num + line] = ''  # Remove processed content make it faster
        temp[0] = temp[0][1:]  # Remove header from sequence label
        output_content.append(temp)
    return output_content  # Return a list with record in the file

# Read in a single line fasta file (2 lines per record)
# FASTA file has to be single line (sequences can not be wrapped into multiple lines)
def read_fasta(input_file):
    fasta = []
    from itertools import zip_longest
    with open(input_file) as f:
        for line1, line2 in zip_longest(*[f]*2):
            fasta.append((line1.strip('\n')[1:], line2.strip('\n')))
    return tuple(fasta)

# Same as read_fasta, but read in FASTQ file (4 lines per record)
def read_fastq(input_file):
    fastq = []
    from itertools import zip_longest
    with open(input_file) as f:
        for line1, line2, line3, line4 in zip_longest(*[f]*4):
            fastq.append((line1.strip('\n')[1:], line2.strip('\n'), line3.strip('\n', line4.strip('\n'))))
    return tuple(fastq)

# Read in a multiple line fasta file, this is much slower than parsing single line FASTA/FASTQ file.
def read_fasta_multiline(input_file, head_symbol='>'):
    input_content = read_file(input_file)
    corrected_content = []
    for line in input_content:
        if line[0] == head_symbol:
            corrected_content.append([line.strip('\n')[1:], ''])
        else:
            corrected_content[-1][1] += line.strip('\n')
    return corrected_content

def read_multiline_fasta(input_file, check_point = 1000):
    fna = []
    i = 0
    with open(input_file, 'r') as f:
        for line in f:
            if line[0] == '>':
                fna.append([line[1:], ''])
                i += 1
                if i >= check_point and i // check_point == 0:
                    print('Loading {0} records ...'.format(i))
            else:
                fna[i-1][1] += line.strip('\n')
    return fna
# Create a Python dictionary for sequencing data (only support FASTA for now)
# By default, label is cut at the first SPACE
# Redundancy is NOT checked at this version.
def create_fna_db(input_content, trim_space ='True', delimiter = ' '):
    db = {}
    for record in input_content:
        if trim_space:
            label = record[0].split(' ')[0]
        else:
            label = record[0]
        db[label] = record[1]
    return db
def write_fna_db(db, db_file):
    import json
    with open(db_file, 'w') as f:
        json.dump(db, f)
    return db_file
def read_fna_db(db_file):
    import json
    with open(db_file, 'r') as f:
        return json.load(f)

def write_json(db, db_file):
    import json
    with open(db_file, 'w') as f:
        json.dump(db, f)
    return db_file

def read_json(db_file):
    import json
    with open(db_file, 'r') as f:
        return json.load(f)

# Create a GFF database so gffutils can read in as database.
def create_gff_db(gff_filename, db_filename):
    import gffutils
    db = gffutils.create_db(gff_filename, db_filename, merge_strategy="merge", force=True)
    print(gff_filename)
    # Print all feature types and count
    print('Feature\tCount\n-------\t-----')
    for item in db.featuretypes():
        print('{0}\t{1}'.format(item, db.count_features_of_type(item)))
    return db

# Read in GFF database created by gffutils
def read_gff_db(db_filename, report = False):
    import gffutils
    db = gffutils.FeatureDB(db_filename)

    if report:
        # Print all feature types and count
        print('Feature\tCount\n-------\t-----')
        for item in db.featuretypes():
            print('{0}\t{1}'.format(item, db.count_features_of_type(item)))
    return db
#%%
# Return reverse compliment of a sequence
# This part is got from Stakoverflow
#(https://stackoverflow.com/questions/19570800/reverse-complement-dna) by corinna
# Works for Python 3
def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhdN', 'TGCAtgcaYRKMyrkmBVDHbvdhN'))[::-1]

# Get sequence by their positions
# start is the starting nucleotide of the gene, but for Python,
# You need to go on character upstream (start-1).
def seq_slice(sequence, start, stop, reverse_compliment = False):
    seq = sequence[start-1:stop]
    if reverse_compliment:
        return revcomp(seq)
    else:
        return seq
#%% Sequence regular expression
# Convert single base to regular expreesion
def convert_IUPAC_to_re(base):
    base = base.upper() # Convert base to upper case
    IUPAC_dict = {'R':['A', 'G'], 'Y':['C', 'T'], 'S':['G', 'C'],\
                  'W':['A', 'T'], 'K':['G', 'T'], 'M':['A', 'C'],\
                  'B':['C', 'G', 'T'], 'D':['A', 'G', 'T'],\
                  'H':['A', 'C', 'T'], 'V':['A', 'C', 'G'],\
                  'N':['A', 'C', 'G', 'T']}
    if base in ['A', 'C', 'G', 'T']:
        return base
    else:
        try:
            return '[' + ''.join(IUPAC_dict[base]) + ']'
        except KeyError:
            print('{0} is not in the IUPAC code.'.format(base))

# Generate the text string for regular expression
def iupac_re(seq):
    seq_se = ''
    for base in seq:
        seq_se += convert_IUPAC_to_re(base)
    return seq_se

# With agiven sequence, if the location is reported on its reverse compliment strand,\
    # return the location of the original strand
# The location should start from 1 instead of 0 for consistence.
def convert_rc_loc(seq_l, rev_start, rev_stop):
    return (seq_l - rev_stop + 1, seq_l - rev_start + 1)

#%% Manipulation on GFF database created with gffutils

#%% Short Oligo search
# This function has several limitations:
    # 1 No gap (since we are dealing with primers, usually gap meaning a bad design)
    # 2 Slow search for error tolerance (Very fast for perfect match with Regular expression,
    # but slower if you want to match with error, like 1 SNP)
    # 3 It will return error if more than one match is found, this is consider that the \
    # primer is not specific, but you should check your sequence to see if you have repeated gene.
    # This function is initially written for rRNA, but should work as a univeral oligo search.

    # It returns a list [start, end, matched sequence], the start position is the real pos on the sequence
    # So it start from 1 (instead of 0)

# A simple search function using regular expression, no mismatch allowed, only search on one strand
# return the matched sequence on the target
# this function is useless for now.
def search(seq, oligo):
    import re
    match = re.findall(iupac_re(oligo), seq, re.I)

    if len(match) > 1:
        return match
    elif len(match) == 0:
        return []
    else:
        return match

# Generate all possible oligo with allowing number of mismatch
# Although it is possible to use mismatch > 1, it is not recommended.
def mismatch_group(oligo, mismatch = 1):
    from itertools import combinations

    if len(oligo) <= mismatch:
        print('Too many mismatches for this oligo, are you crazy?')
    elif mismatch == 0:
        return [oligo]

    pos = [x for x, y in enumerate(oligo)] # Put position of bases into a list
    mismatch_oligos = []
    for item in combinations(pos, mismatch):
        oligo_alternative = ''
        for x, y in enumerate(oligo):
            if x in item:
                oligo_alternative += 'N'
            else:
                oligo_alternative += y
        mismatch_oligos.append(oligo_alternative)
    return mismatch_oligos

# Search a DNA sequence with oligos. Mismatch is allowed with any number.
# Larger than 2 mismatches will be super slow.
# Do NOT support indels for now.
# All possible matches will be reported, lower case indicated the match is on\
    #reverse strand.
# 2 indicate the search is done on both strands, 1 indicate the search is on\
    #the given strand (can be plus or minus).

# Always need re.I since orginal Genome maybe marked with lower case bases.
def search_oligo(seq, oligo, mismatch = 0, rev_strand = True):
    import re
    match_list =[]
    oligo_alternatives = [oligo]
    oligo_alternatives += mismatch_group(oligo, mismatch = mismatch)
    #print(oligo_alternatives)
    for item in oligo_alternatives:
        match = re.findall(iupac_re(item), seq, re.I)
        match_list += [i.upper() for i in match]
    # Allowing mismatch will result in finding same sequencing using different alternatives
    # Solve this using set
    if mismatch > 0 and len(match_list) > 1:
        temp_list = set(match_list)
        if len(temp_list) < len(match_list): # Found redundnat sequences
            match_list = []
            for item in temp_list:
                match_list += (re.findall(item, seq, re.I))

    if not rev_strand:
        return [match_list, 1]
    # Search again on reverse compliment strand (same procedure).
    else:
        seq_rev = revcomp(seq)
        match_list_rev = []
        for item in oligo_alternatives:
            match = re.findall(iupac_re(item), seq_rev, re.I)
            match_list_rev += [i.lower() for i in match]
        # If mismatch is allower, matching can be redudant (from different primer alternates)
        # We need to remove the redundancy before proceed.
        if mismatch > 0 and len(match_list_rev) > 1:
            temp_list = set(match_list_rev)
            if len(temp_list) < len(match_list_rev): # Found redundnat sequences
                match_list_rev = []
                # After dereplicate the searching result, we search again use the non-degenrate\
                # sequence (which is much faster)
                for item in temp_list:
                    match_list_rev += (re.findall(item, seq_rev, re.I))
        return [[i.upper() for i in match_list] + [i.lower() for i in match_list_rev], 2]

# Return the location of a primer in a given sequences
# The search will be on both strand,\
# reported location is on the given strand
# Since we are using the result of search_oligo, mismatch has already been dealed with
# Starting locatoin is the real location on sequence, not for Python slicer.
# All locations are reported as on the original strand.
def search_location(seq, oligo):
    import re
    #seq_rev = revcomp(seq)
    loc_set = []
    for loc in re.finditer(oligo, seq, re.I):
        loc_set.append((loc.start() + 1, loc.end())) # Return the location as a tuple
    return loc_set

# A combined function from primer search to primer loci.
def search_oligo_loc(seq, oligo, mismatch = 0, rev_strand = True):
    # Get all candidate match on target
    m = search_oligo(seq, oligo, mismatch = mismatch, rev_strand = rev_strand)
    loc = {}
    for item in m[0]:
        if item.isupper(): # target is on original strand
            loc[item] = search_location(seq, item)
        elif item.islower(): # target is on the other strand
            loc[item] = [convert_rc_loc(len(seq), i[0], i[1]) for i in search_location(revcomp(seq), item)]
        else:
            print(item)
            print('You have mixed UPPER and lower case letter in your data, please check.')
    return loc
#%% PCR
# With given primers, return the amplicon sequence or a report if no match or multiple match.
def pcr(seq = '', pf= '', pr = '', mismatch_pf = 0, mismatch_pr = 0):
    f_match = search_oligo(seq, pf, mismatch = mismatch_pf, rev_strand=False)
    r_match = search_oligo(revcomp(seq), pr, mismatch = mismatch_pr, rev_strand=False)
    print([f_match, r_match])

    if len(f_match[0]) == 1 and len(r_match[0]) == 1: # Both primers have unique match
        start = f_match[0][0]
        stop = len(seq) - r_match[0][0] + 1
        amplicon_seq = seq_slice(seq, start, stop)
    return amplicon_seq

#%% Simulate DNA fragmentation
# If we assume DNA fragmentatin following normla distribution
# That is: an average fragment length, with the variation of lengtg
def fragment_dist(seq_length = 5000, mean = 200, std = 1000):
    import numpy as np
    fragment_list = []
    while seq_length > 0:
        random_fragment = int(np.random.normal(mean, std))
        if random_fragment > 0:
            if random_fragment <= seq_length:
                fragment_list.append(int(random_fragment))
                seq_length -= int(random_fragment)
            if random_fragment > seq_length:
                fragment_list.append(seq_length)
                seq_length -= int(random_fragment)
    return fragment_list

# Generate the DNA fragment from a Given sequence using the fragment length list
# Fragment shorter or longer than expected length were filter out
# Fragment were picked from the 5' of the given sequence
def generate_fragment(seq, mean = 200, std = 500, filter_low = 50, filter_high = 1000):
    pass
    # return a FASTA list


#%% Simulate DNA mutation
# Generate SNP on the given DNA sequences based on similarities
# If similarity is 0.97, then it will generate a sequence with 3% difference
# SNP loci are draw randomly without replacement, so there is chance for two or three continuous SNPs.
def mutate_seq(seq, similarity):
    import random
    snp_dict = {'A':['T','C','G'], 'T':['A','C','G'],\
                'C':['A','T','G'], 'G':['A','T','C']}
    seq = seq.upper()
    seq = [i for i in seq]
    if similarity > 1:
        print('Similarity need to be between 0 and 1.')
        return None
    else:
        #print(len(seq) * (1-similarity))
        pos = range(0, len(seq))
        pos_sample = random.sample(pos, int(len(seq) * (1 - similarity)))
        for item in pos_sample:
            seq[item] = random.sample(snp_dict[seq[item]], 1)[0]
        seq = ''.join(seq)
        return seq


#%% Kmer handlings
# Several slow kmer functions for testing
# May need to consider using khmer on real pipelines


# Return the kmer list of a given sequences
# It will return all kmer in both strand WITHOUT dereplication
# So the kmer will be REDUNDANT in the output
def kmer(seq, size):
    seq = seq.upper()
    seq_length = len(seq)
    kmer = []
    for i in range(seq_length - size + 1):
        kmer.append(seq[i:i+size])
    seq_rev = revcomp(seq)
    for i in range(seq_length - size + 1):
        kmer.append(seq_rev[i:i+size])
    return kmer


# Return a kmer table with abundance on the given kmer length for a set of sequences
# The input sequences need to be in the style as the otuput of read_seqs() or read_fasta(),
# In which sequences are store in tuple as (label, sequecne), (label, sequence), ...)
def kmer_count(seqs, size):
    kmer_list = []
    for item in seqs:
        kmer_list += kmer(item[1], size)
    kmer_table = {}
    for item in kmer_list:
        kmer_table[item] = kmer_table.get(item, 0) + 1 # one line code for adding new keys
#        try:
#            kmer_table[item] += 1
#        except KeyError:
#            kmer_table[item] = 1
    return kmer_table


# Return the kmer table by compare two kmer table from kmer_table().
# In order to speed up, we use the kmer_count as input, instead of raw sequences
def kmer_pair(count1, count2, size):
    kmers = set(list(count1.keys()) + list(count2.keys()))
    kmer_table = []
    for item in kmers:
        kmer_table.append((item, count1.get(item, 0), count2.get(item, 0)))
    kmer_table.sort(key=lambda x:x[0])
    return kmer_table


# Return the kmer distance, Jaccard distance, NOT account for kmer abundance
def kmer_jaccard(kmer_table):
    total_size = len(kmer_table)
    shared_size = 0
    for line in kmer_table:
        if line[1] > 0 and line[2] > 0: # Current kmer is present in both sets
            shared_size += 1
    return shared_size / total_size

def kmer_mash_distance(jaccard, size):
    import numpy as np
    jac = jaccard
    if jac == 0:
        return 1
    else:
        mash_distance = (-1/size) * np.log((2 * jac)/(1 + jac))
        return mash_distance


# Return the kmer distance, accounted for kmer abundance, usng Euclidean distance
# May need to standardize before calculating
def kmer_euclidean(set1, set2, size):
    kmer_table = kmer_pair(set1, set2, size)
    distance = sum([(i[1] - i[2]) ** 2 for i in kmer_table]) ** 0.5
    return distance