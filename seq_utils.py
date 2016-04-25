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

#####################
# GENERAL FUNCTIONS #
#####################

def quantile_norm(expression_list):
    """
    Given a list,
    returns that list normalized by quantiles after performing a log2 normalization
    """

def merge_expression_experiments(directory, quant=True):
    """
    Given a directory with n pile files
    generates a matrix file with

    position exp1 exp2 exp3 .... exp n-1 exp n

    The expression is passed to log2 and normalized by quantiles (this last can be prevented by arguments)
    """



#####################
#      CLASSES      #
#####################


#####################
#     EXECUTION     #
#####################

if __name__ == "__main__":

