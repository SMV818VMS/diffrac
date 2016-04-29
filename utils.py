#/usr/bin/env python

#############################################################
#
# utils.py
#
# Author : Miravet-Verde, Samuel
# Written : 04/25/2016
# Last updated : 04/25/2016
#
# General utils functions
#
#############################################################

#####################
#   PACKAGES LOAD   #
#####################

import sys, os
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import seaborn
from scipy.signal import savgol_filter

#####################
# GENERAL FUNCTIONS #
#####################

def return_column(filename, index, header=True):
    """
    Given a file with columns
    returns a list with the numbers appearing in tha column
    """

    fi = open(filename, 'r')
    column = []

    for line in fi:
        if header:
            header = False
        else:
            column.append(float(line.strip().split()[index]))

    fi.close()
    return column


def regression(x, y, limit, title='decay fit', location=True, soft=True):
    """
    Performs a basic regression fitting giving two ARRAY structures x and y

    The location of the plot can be determined

    Returns slope and value
    """


    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

    if location:
        # Define the regression line
        ablineValues = []
        for i in x[:limit]:
            ablineValues.append(slope*i+intercept)

        # Plot points and fitting line
        plt.plot(x, savgol_filter(y, 21,4), c='k', label='expression', linewidth=1.8, alpha=0.7)
        plt.plot(x[:limit], ablineValues, c='r', label='linear fit', linewidth=1.5, alpha=0.8)
        plt.title(title)

        # Add annotation
        plt.text(x[0], max(y), 'r2='+str(r_value**2))
        plt.text(x[0], max(y)-10, 'slope='+str(slope))

        plt.savefig('results_finddrops/'+title+'.png')
        plt.close()

    return(slope, r_value**2)



#####################
#      CLASSES      #
#####################



#####################
#     EXECUTION     #
#####################

if __name__ == '__main__':
    regression([1,2,3,4,5,6,7,8,9,10],[10,10,10,10,10,10,10,10,10,0], title='sharp')
    regression([1,2,3,4,5,6,7,8,9,10],[10,6,3,2,1,0.5,0,0,0,0], title='decay1')
    regression([1,2,3,4,5,6,7,8,9,10],[10,9.8,9.5,9,8,5,2,0.5,0,0], title='decay2')
