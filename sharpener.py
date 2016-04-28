#/usr/bin/env python

#############################################################
#
# sharperner.py
#
# Author : Miravet-Verde, Samuel
# Written : 04/18/2016
# Last updated : 04/228/2016
#
# A different approach to detect sharp changes in expression
#
#############################################################

#####################
#   PACKAGES LOAD   #
#####################

import sys, os
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import utils as u
import seq_utils as sequ
import seaborn

#####################
# GENERAL FUNCTIONS #
#####################


def drop_stats():
    """
    Function required to know if a drop in expression can be considered a sharp
    change or not.

    This algorithm looks in all the annotation, to define which is a significant change
    to consider a drop between contiguous positions
    """

    #1. Extract the expression profile for all the available annotations
    expression  = sequ.annotated_expression(log=False)

    #2. First derivative
    dexpression = np.diff(expression)
    dexpression_neg = [abs(x) for x in dexpression if x < 0]
    dexpression_pos = [x for x in dexpression if x > 0]

    #3. Basic stats
    mean_neg, mean_pos = np.mean(dexpression_neg), np.mean(dexpression_pos)
    stdv_neg, stdv_pos = np.std(dexpression_neg), np.std(dexpression_pos)

    std2neg = mean_neg + (stdv_neg*2)
    std2pos = mean_pos + (stdv_pos*2)

    # Pass to log2
    dexpression_neg = [np.log2(x) for x in dexpression_neg]
    dexpression_pos = [np.log2(x) for x in dexpression_pos]

    #4. Plot the distribution
    bins = 50
    plt.hist(dexpression_neg, bins, alpha=0.5, label='negative')
    plt.hist(dexpression_pos, bins, alpha=0.5, label='positive')
    plt.xlabel('log2(change)', fontsize=18)
    plt.ylabel('frequency', fontsize=18)
    plt.title('log2 change between contiguous base positions', fontsize=18)
    plt.legend(loc='upper right', fontsize='large')

    plt.axvline(x=np.log2(mean_neg), c='b', linestyle='--', label=r'$\mu\ negative$', linewidth=2, alpha=0.5)
    plt.axvline(x=np.log2(mean_pos), c='g', linestyle='--', label=r'$\mu\ positive$', linewidth=2, alpha=0.5)
    plt.axvline(x=np.log2(stdv_neg), c='b', linestyle=':', label=r'$\mu+2\sigma$', linewidth=2)
    plt.axvline(x=np.log2(stdv_pos), c='g', linestyle=':', label=r'$\mu+2\sigma$', linewidth=2)

    plt.savefig('./figures/sharpener_hist.png')
    plt.close()


def sharpener(window):
    """
    Given a window with expression values
    Discriminates if it can be a sharp drop or not by means the
    first derivative
    ________             ________ ________
            |____ ---->          |         ----> Sharp
    """

    
#####################
#      CLASSES      #
#####################


#####################
#     EXECUTION     #
#####################

if __name__ == "__main__":
    drop_stats()
