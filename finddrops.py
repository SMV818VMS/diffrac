#/usr/bin/env python

#############################################################
#
# finddrops.py
#
# Author : Miravet-Verde, Samuel
# Written : 04/18/2016
# Last updated : 04/20/2016
#
# Test script to detect drops in the expression/toyset files
#
#############################################################

#####################
#   PACKAGES LOAD   #
#####################

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import seaborn

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


def decay_score(your_list, index):
    """
    Define the decay score for a list of numbers

    The decay score include to numbers:
    [stdsc, maxsc, dropsc, decaysc]

    stdsc   = standard deviation for the expression in that window
    maxsc   = minimum number in the first derivative function, i.e. maximum drop
    dropsc  = difference between last expression value and first 0
    decaysc = the decay start that marks the first nucleotide after which the expression profiles starts to decay with the factor

    This two values allow to differentiate between a sharp and a decay termination:
    --> high stdsc, low maxsc and maxsc != dropsc ==> decay
        The change in expression decays gradually but with no big changes, this hardens the match between maximum change (minimum derivative)
        and the drop in the index position
    --> low stdsc, high maxsc and maxsc == dropsc ==> sharp
        The change in expression is abrupt and with a big difference, this makes easier the match between minimum derivative and the index position.

    Additionally, if the termination is sharp, the decaysc will match the last expression base, thing that will not occur if the termination
    is in decay.
    """

    a = np.std(your_list[:index+1])                    # +1 include last expression value
    first_derivative = np.diff(your_list[:index+2])     # +2 include first 0
    b = min(first_derivative)
    c = list(first_derivative)[index]

    # Define the position where the decay starts:
    # 1. Find the last positive value in the first derivative:
    positive_indexes = [i for i, e in enumerate(first_derivative) if e > 0]
    # 2. The maximum will correspond to the last position, we have to sum one to point to the
    # first nucleotide after which the expression only decreases
    decay_start = max(positive_indexes)+1

    return([a, b, c, decay_start])


def plot_drop(start, your_list, main_title, last_expression, decay_start):
    """
    Plot the drop
        start = defines the first position for the x label
        your_list = the list of point to plot
        main_title = title for the plot
        last_expression = defines the position where the non expression starts
        decay_start =  defines the position with the last positive 1st derivative (after this, the expression drops)
    """

    # Define the data we want to plot
    x = list(range(start, start+len(your_list)))
    y = your_list

    # Define the figure, the labels and plot it
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax1.set_title(main_title)
    ax1.set_xlabel('genome position (pb)')
    ax1.set_ylabel('expression counts')

    ax1.plot(x, y, c='k', label='expression', linewidth=2)
    ax1.axvline(x=last_expression, c='b', label='last_expression', linewidth=2)
    ax1.axvline(x=decay_start, c='r', label='decay_start', linewidth=2)

    leg = ax1.legend(fontsize='large')

    plt.show()


def find_drops(annotation_file, expression_file, expression_index, expression_threshold=0.0, expression_determinant=4, decay_window=100, header_ann=True, header_exp=True):
    """
    Given a pile up file,
    returns the positions with potential termination signals

    Parameters:
    annotation_file:        file including the annotation of the genome (TSS and TTS), used to define the window of work
    expression_file:        non fractionates RNA Seq file
    expression_index:       index of the column with the expression values in the expression file
    expression_threshold:   value in expression to consider a position non expressed
    expression_determinant: factor that divides the intergenic region to delimit the window
    decay_window:           window previous to the non expressed window considered
    header_ann:             boolean, define if the annotation file has headers or not
    header_exp:             boolean, define if the expression file has headers or not
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
    c = 1

    while i < len(expression)-sliding_window:
        # Define the window of work
        current_window = expression[i:i+sliding_window]

        # Only analyze the window if the expression drops below the threshold in the no_exp_window after a value with expression:
        if current_window[:decay_window+1].count(expression_threshold) == 0.0 and current_window[decay_window-1] != 0 and np.mean(current_window[decay_window+1:]) <= expression_threshold:
            stdsc, maxsc, dropsc, decaysc = decay_score(current_window, decay_window)
            identifier = 'SIGN'+str(c)
            last_expression = i+decay_window+1
            decay_start = i+decaysc
            results[identifier] = [i, i+sliding_window, decay_start, last_expression, stdsc, maxsc, dropsc]
            c += 1

            # Plot the drop
            plot_drop(i, current_window, identifier, last_expression, decay_start)
        i+=1

    print(results)


#####################
#      CLASSES      #
#####################


#####################
#     EXECUTION     #
#####################

if __name__ == "__main__":
    find_drops(annotation_file='./datasets/toyset_annotations.txt', expression_file='./datasets/toyset.txt', expression_index=1, header_exp=False)
