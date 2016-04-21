#/usr/bin/envimport sys, os

#############################################################
#
# findmotifs.py
#
# Author : Miravet-Verde, Samuel
# Written : 04/20/2016
# Last updated : 04/21/2016
#
# Test to detect motifs from sequences given by finddrops
#
#############################################################

#####################
#   PACKAGES LOAD   #
#####################

import sys, os
import itertools
import numpy as np
import glob
from scipy.cluster.vq import kmeans, vq
from finddrops import find_drops
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

#####################
# GENERAL FUNCTIONS #
#####################

def maxgap(data, max_gap):
    """
    Arrange data into groups where successive elements
    differ by no more than *maxgap*

    The input has to be a list of list with the structure:
    [['id1', distance between start decay and last expression],['id2',dist2],...,['idn', distn]]

    The output will be a list of lists with the identifiers clustered together by the distances if
    the difference between them is less than the maxgap
    [[id1, id2],[id3,..]]

    Example:
    in: [['id1',1], ['id2',-1], ['id3',2], ['id4',80], ['id5',81], ['id3',82]]
    out: [['id1','id2','id3'], ['id4','id5','id3']]
    """

    # Sort the list by the second value (distances)
    data.sort(key=lambda x: x[1])

    # Separate lists
    identifiers = [x[0] for x in data]
    distances   = [x[1] for x in data]

    # Cluster the distances by max gap and return the identifiers groupes
    groups = []
    for k, g in itertools.groupby(distances, key=lambda n: n//max_gap):
        i = len(list(g))
        groups.append(identifiers[:i])
        identifiers = identifiers[i:]
    return groups


def group_by(data, method, max_gap=10, k=2):
    """
    Define the number of groups in your dataset and labels them
    Possible methods:
    - maxgap
    - kmeans

    ----

    max_gap: maximum distance expected between features in maxgap method
    k: is the number of groups expected, used by some of the methods
    """

    if method == 'maxgap':
        return maxgap(data, max_gap)


def subset_dictionary(bigdict, wanted_keys):
    """
    Given a dictionary (bigdict) including several keys and a list of keys
    in list format in wanted_keys

    Return a subsetted dictionary from bigdict with only wanted keys
    """

    return dict((k, bigdict[k]) for k in wanted_keys if k in bigdict)


def process_sequences(dictionary, genome, outFile):
    """
    Given a dictionary with {identifier:[features]}
    where feature[0] is the start index of a sequence and feature[1] is the end

    Creates a multifasta file to run a meme analysis
    """

    # The goal is extract several subsequences from a genome given
    # we load the genome

    with open(outFile, 'w') as f:
        for seq_record in SeqIO.parse(genome, 'fasta'):
            for k, v in dictionary.iteritems():
                start = v[0]
                end   = v[1]
                f.write('>'+k+'_'+str(seq_record.id)+'\n')
                f.write(str(seq_record.seq[start:end])+'\n')


def pymeme(directory):
    """
    Given a directory,
    List all the fasta files in it and runs a MEME for each file

    By default, the program creates a directory with the name of the fasta where
    all the results will be located
    """

    print('finding files...')

    files = glob.glob(directory+"*.fasta")

    print(str(len(files))+' files found')

    for fil in files:
        outname = fil.replace(directory,'').replace('.fasta','')
        print('meme of '+outname)
        command = 'meme '+fil+' -dna -oc ./results_meme/'+outname+'/ -nostatus -time 18000 -maxsize 60000 -mod zoops -nmotifs 5 -minw 6 -maxw 100'
        os.system(command)


#####################
#      CLASSES      #
#####################



#####################
#     EXECUTION     #
#####################

if __name__ == "__main__":

    # Generate the list of candidates
    candidates = find_drops(annotation_file='./datasets/toyset_annotations.txt', expression_file='./datasets/toyset.txt', expression_index=1, header_exp=False)

    # This is a dictionary, to group we have to select a feature explaining more differences, in this dictionary we have different features:
    # 'SIGN7': [6986, 7211, 7001, 7087, 48.566111083270691, -68.0, -4.6440034020500001]
    # ident : [start, end,decay_p, 0exp, std expression, max_change, diff exp btw last and 0exp]

    # The most specific value for classes is the difference between decay position and 0exp
    # Reformat the dictionary

    data = []
    for key, values in candidates.iteritems():
        distance = values[3]-values[2]
        data.append([key, distance])

    # Group by
    groups = group_by(data, 'maxgap')

    # Extract subsequences and write them in separate files:
    sharp_signals = subset_dictionary(candidates, groups[0])
    decay_signals = subset_dictionary(candidates, groups[1])

    process_sequences(sharp_signals, './datasets/toy_genome.fasta', './results_meme/input_sequences/sharp_signals.fasta')
    process_sequences(decay_signals, './datasets/toy_genome.fasta', './results_meme/input_sequences/decay_signals.fasta')

    # Run pymeme
    pymeme('./results_meme/input_sequences/')
