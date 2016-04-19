#/usr/bin/env python

#############################################################
#
# finddrops.py
#
# Author : Miravet-Verde, Samuel
# Written : 04/18/2016
# Last updated : 04/19/2016
#
# Test script to detect drops in the expression/toyset files
#
#############################################################

#####################
#   PACKAGES LOAD   #
#####################

import sys, os
import numpy as np

#####################
# GENERAL FUNCTIONS #
#####################

def intergenic_mean(filename, header=True):
    """
    Given an annotation dataset,
    Returns the mean distance between genic regions

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
    previous_end = None
    number_ann = 0
    sum_distan = 0

    with open(filename, 'r') as fi:
        for line in fi:
            if header == True:
                header = False
            else:
                gene, current_start, current_end = line.strip().split()

                # Compute mean without considering the overlapping genes
                if previous_end and previous_end < int(current_start):
                    number_ann += 1
                    sum_distan += (int(current_start)-int(previous_end))

                previous_end = int(current_end)

    return(sum_distan/float(number_ann))


def return_column(filename, index, header=True):
    """
    Given a file with columns
    Returns a list with the column determined by index (starting in 0)
    """

    fi = open(filename, 'r')
    column = [float(line.strip().split()[index]) for line in fi]
    fi.close()
    return(column)


def decay_score(your_list):
    """
    Define the decay score for a list of numbers

    The decay score include to numbers:
    [a, b]

    a = mean of changes for the first derivative function
    b = minimum number in the first derivative function

    This two values allow to differentiate between a sharp and a
    decay termination:

    --> high a and low b = decay
        The change in expression decays gradually but with no big changes
    --> low a and high b = sharp
        The change in expression is abrupt and with a big difference
    """

    first_derivative = np.diff(your_list)
    a = np.mean(first_derivative)
    b = min(first_derivative)

    return([a, b])


def find_drops(annotation_file, expression_file, expression_index, expression_threshold=0.0, expression_determinant=4, decay_window=100, header_ann=True, header_exp=True):
    """
    Given a pile up file,
    returns the positions with potential termination signals
    """

    # 1. First we have to define the mean distance between annotations to set up a correct
    # window size to run the algorithm

    intergenic_distance_mean = intergenic_mean(annotation_file)


    # 2. With this value we have the expected size of regions with expression equal 0.
    # The algorithm runs windows along the genome trying to detect these regions falling to 0 and looking a
    # defined number of bases before and computes a decay factor.
    #      - This decay factor is computed as the number of bases we need to pass from the maximum expression
    #      value in the window until we arrive to 1/4 intergenic distance number of zeros.

    # Define the windows to process the expression profile
    no_exp_window  = int(round(intergenic_distance_mean/expression_determinant))
    sliding_window = int(round(decay_window + no_exp_window))

    print(no_exp_window, decay_window, sliding_window)

    # Load the expression profile and start the sliding window process
    # While processing, append the results to the dictionary of results
    results = {}
    expression = return_column(expression_file, expression_index)

    i = 0
    c = 0
    while i < len(expression)-sliding_window:
        # Define the window of work
        current_window = expression[i:i+sliding_window]

        # Only analyze the window if the expression drops below the threshold in the no_exp_window after a value with expression:
        if current_window[:decay_window+1].count(expression_threshold) == 0.0 and current_window[decay_window-1] != 0 and np.mean(current_window[decay_window+1:]) <= expression_threshold:
            print(current_window)
            c += 1
        else:
            pass

        i+=1

    print(c)

#####################
#      CLASSES      #
#####################


#####################
#     EXECUTION     #
#####################

if __name__ == "__main__":
    find_drops(annotation_file='./datasets/toyset_annotations.txt', expression_file='./datasets/toyset.txt', expression_index=1, header_exp=False)
