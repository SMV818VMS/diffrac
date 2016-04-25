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


#####################
#      CLASSES      #
#####################



#####################
#     EXECUTION     #
#####################

