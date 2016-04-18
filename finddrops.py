#/usr/bin/env python

#############################################################
#
# finddrops.py
#
# Author : Miravet-Verde, Samuel
# Written : 04/18/2016
# Last updated : 04/18/2016
#
# Test script to detect drops in the expression/toyset files
#
#############################################################

#####################
#   PACKAGES LOAD   #
#####################

import sys, os
from itertools import tee
from itertools import izip

#####################
# GENERAL FUNCTIONS #
#####################

def pairwise(iterable):
    """ s -> (s0,s1),(s1,s2),(s2,s3)...(sN-1,sN)"""
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)





def intergenic_mean(filename, header=True):
    """
    Given an annotation dataset,
    Rerturns the mean distance between genic regions

    Ex:
    < file
    gene TSS TTS
    1    1   2
    2    5   10
    3    16  19


    The program runs:
    (5-2)+(16-10) = 3+6 = 9
    9/2 = 4.5
    """

    annotations = []

    with open(filename, 'r') as fi:
        for line in fi:
            if header == True:
                header = False
            else:
                l.append()






#####################
#      CLASSES      #
#####################


#####################
#     EXECUTION     #
#####################


# 1. First we have to define the mean distance between annotations to set up a correct
# window size to run the algorithm


