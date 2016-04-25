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
import utils as u
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

def decay_score(your_list, index, wooble=0.0):
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

    Wooble is dependant of the noise in the expression profile. It is considered to determine the decaysc as a 1 or a 2 in the decay tail
    can be considered as the first with positive expression but it comes from the oscllation of the profile.
    """

    stdsc = np.std(your_list[:index+1])                    # +1 include last expression value
    first_derivative = np.diff(your_list[:index+2])     # +2 include first 0
    maxsc = min(first_derivative)
    dropsc = list(first_derivative)[index]

    # This is more tricky...
    # Define the position where the decay starts:
    # 1. Find the last positive value in the first derivative:
    # 2. The maximum will correspond to the last position, we have to sum one to point to the
    # first nucleotide after which the expression only decreases

    positive_indexes = [i for i, e in enumerate(first_derivative) if (e-wooble) > 0.0]
    decay_start = max(positive_indexes)+1

    return([stdsc, maxsc, dropsc, decay_start])


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

    plt.savefig('./results_finddrops/'+main_title+'.png')

    plt.close()


def find_drops(annotation_file, expression_file, expression_index, expression_threshold=0.0, decay_variability = 0.0, expression_determinant=4, decay_window=200, header_ann=True, header_exp=True):
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

    # Load the expression profile and start the sliding window process
    # While processing, append the results to the dictionary of results
    results = {}
    expression = u.return_column(expression_file, expression_index)

    i = 0
    c = 1

    while i < len(expression)-sliding_window:
        # Define the window of work
        current_window = expression[i:i+sliding_window]

        # Only analyze the window if the expression drops below the threshold in the no_exp_window after a value with expression:
        if current_window[:decay_window+1].count(expression_threshold) == 0.0 and current_window[decay_window-1] != 0 and np.mean(current_window[decay_window+1:]) <= expression_threshold:
            stdsc, maxsc, dropsc, decaysc = decay_score(current_window, decay_window, decay_variability)
            identifier = 'SIGN'+str(c)
            last_expression = i+decay_window+1
            decay_start = i+decaysc+1
            results[identifier] = [i, i+sliding_window, decay_start, last_expression, stdsc, maxsc, dropsc]
            c += 1

            # Plot the drop
            plot_drop(i, current_window, identifier, last_expression, decay_start)
        i+=1

    # Write the file with the results:
    # iterate by keys in sorted order
    fo = open('./results_finddrops/drop_signals.txt','w')
    fo.write('id\tstart\tend\tdcy_st\tlast_exp\tstdsc\tmaxsc\tdropsc\n')
    for k in sorted(results.keys()):
        value = '\t'.join([str(i) for i in results[k]])
        fo.write(k+'\t'+value+'\n')

    fo.close()
    return results

#####################
#      CLASSES      #
#####################


#####################
#     EXECUTION     #
#####################

if __name__ == "__main__":
#    find_drops(annotation_file='./datasets/toyset_annotations.txt', expression_file='./datasets/toyset.txt', expression_index=1, header_exp=False)
    find_drops(annotation_file='../mycorepo/plusTSSTTS.csv', expression_file='./datasets/dsspilesmpn.txt', expression_index=2, decay_variability=10.0, expression_determinant=12)
