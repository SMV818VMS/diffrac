#!/usr/bin/env python

from random import choice
from random import randint as r
from math import exp as e
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

###########################################################
# EXPRESSION FILE
###########################################################

# The idea is to create a toy dataset including the two types of termination
# we can see in the RNA seq profiles to test two different things:

# 1. The size of the windows required to detect the signals
# 2. The whole approach

max_exp = 200
genome_size = 9000


# We need a decay list of numbers:
decay = []
c = 1/e(2)
x = max_exp

while x >= 0:
    decay.append(x)
    c = c*1.05
    x = x-c


# Arrive to 0
decay = decay[1:]+[0]

# We will put the information in the toyset file
fo = open('./datasets/toyset.txt','w')

# Define the properties of the population and the position for 2 types of
# termination, decay and sharp. We define the list of positions to be able to test they work:

sharp_pos = [2000, 4000, 6000, 8000]
decay_pos = [1000, 3000, 5000, 7000]
window = 100
no_exp = 500

expression = [r(max_exp/2, max_exp) for x in range(genome_size)]

for i in sharp_pos:
    expression[i:i+no_exp] = [0]*no_exp

for i in decay_pos:
    expression[i:i+len(decay)] = decay
    expression[i+len(decay):i+len(decay)+no_exp] = [0]*no_exp


# Create expression file:
for i in range(1, genome_size+1):
    fo.write(str(i)+'\t'+str(expression[i-1])+'\n')

fo.close()

###########################################################
# GENOME FILE
###########################################################

# Create a random genome

genome = ''.join(choice(['A','C','G','T']) for _ in range(genome_size))

# Add motifs to be detected

# Define the motifs:
sharp_motif = 'CCAATCGCTAAATCGCTGTTGGTATTTTTATT'
decay_motif = 'GATAGATAGATAGATAGATAGATAGATAGATA'

# Add the motifs:

def add_motif(genome, motif, index):
    """
    Replace from the index to index+len(motif) in the genome
    with the motif selected
    """

    try:
        # If this can be done is because index is an integer
        index = index-len(motif)
        return genome[:index]+motif+genome[index+len(motif):]
    except:
        # If we arrive here is becuse index is a list
        for i in index:
            i = i-len(motif)
            genome=genome[:i]+motif+genome[i+len(motif):]
        return genome

genome = add_motif(genome, sharp_motif, sharp_pos)
genome = add_motif(genome, decay_motif, decay_pos)

# Create a biopython object to write the sequence in fasta format
handle = open('./datasets/toy_genome.fasta','w')
seq = SeqRecord(Seq(genome, generic_dna), id='ToyGenome')
SeqIO.write(seq, handle, "fasta")
handle.close()
