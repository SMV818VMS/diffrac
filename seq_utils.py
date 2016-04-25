#/usr/bin/env python

#############################################################
#
# seq_utils.py
#
# Author : Miravet-Verde, Samuel
# Written : 04/25/2016
# Last updated : 04/25/2016
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


def kmeans(data, k):
    """
    Arrange data into k groups based on a kmeans approach

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
    print(group_by(data, 'maxgap'))
    print(group_by(data, 'kmeans'))

