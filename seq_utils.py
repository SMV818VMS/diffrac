#/usr/bin/env python

#############################################################
#
# seq_utils.py
#
# Author : Miravet-Verde, Samuel
# Written : 04/25/2016
# Last updated : 04/27/2016
#
# Test to detect motifs from sequences given by finddrops
#
#############################################################

#####################
#   PACKAGES LOAD   #
#####################

import sys, os
import glob
import re
import utils as u
import pandas as pd
import numpy as np


#####################
# GENERAL FUNCTIONS #
#####################

# BASIC SEQUENCE FUNCTIONS

def change_base(i, base=816394):
    """
    Basic function to pass from + to - strand

    1   --> 816394
    300 --> 816094
    """

    return (base+1)-i


# TRANSCRIPTOMIC DATA LOADING AND INTEGRATION

def pseudo_log2(i):
    return np.log2(i+1)


def load_expression(fil, sep='\t', header=0, log=True):
    """
    Given a expression file with two columns, one per strand
    Returns a dataframe of it

    This function also:
    - log2
    - corrects the assignment of the strand considering the rRNA expression
    """

    # Load
    df = pd.read_csv(fil, sep=sep, header=header)

    # Position column to rownames and remove it (also removes a unnamed column)
    df.set_index(df[df.columns[0]], inplace=True)
    try:
        df.drop(df.columns[[3]], axis=1, inplace=True)
    except:
        pass
    df.drop(df.columns[[0]], axis=1, inplace=True)

    # Extract a identifier
    # sample_name = fil.rsplit('/')[-1].rsplit('_',1)[:-1][0]
    sample_name = fil.rsplit('/')[-1].replace('.pile', '')

    # Check rRNA expression and interchange columns if necessary
    # Independently of plus/minus, add a identifier to the colname
    if np.mean(df.iloc[120000:123000,0]) < np.mean(df.iloc[120000:123000,1]):
        # The greater mean is the plus strand
        df.columns = ['minus_'+sample_name, 'plus_'+sample_name]
    else:
        df.columns = ['plus_'+sample_name, 'minus_'+sample_name]

    # Log2 transformation
    if log:
        df = df.apply(pseudo_log2)

    return df


def normalization(df, positive=True):
    """
    Given a df,
    returns the df normalized
    """

    if positive:
        df_norm = (df - df.min()) / (df.max() - df.min())
    else:
        df_norm = (df - df.mean()) / (df.max() - df.min())
    return df_norm


def merge_expression_experiments(directory, norm=True):
    """
    Given a directory with n pile files
    generates a matrix file with

    position exp1 exp2 exp3 .... exp n-1 exp n

    The expression is passed to log2 and normalized by quantiles (this last can be prevented by arguments)
    """

    # List all the files in a directory if they are piles
    files = glob.glob(directory+'*.pile')

    # Load the data iteratively and merge all of them

    first = True

    for fil in files:
        if first:
            expression = load_expression(fil)
            first = False
        else:
            expression = expression.append(load_expression(fil))

    # Normalize
    if norm:
        expression = normalization(expression)

    return expression


# GENOME EXPRESSION INFORMATION


def annotated_expression(log=True):

    #1. Load data and prepare a dictionary
    files = ['/home/smiravet/crg/mycorepo/minusTSSTTS.csv', '/home/smiravet/crg/mycorepo/plusTSSTTS.csv']
    expression = load_expression('./datasets/dsspilesmpn.txt', log=log)
    expressed = []

    #2. Assign mean and standard deviation to all the annotations
    for i in [0,1]:
        with open(files[i], 'r') as fi:
            for line in fi:
                genename, start, end = line.strip().split()
                start, end = int(start), int(end)
                if re.match('^[A-Z]{3}[0-9]{3}$', genename):
                    expressed += list(expression.iloc[start:end, i])


    return expressed


def no_expression_guesser(percentage=5):
    """
    Define which is the expected value to consider the expression is dropped

    To compute it, this takes the mean expression value for codifying regions
    and returns the mean + 1Standard Deviation

    If percentage is set to 0, it takes into consideration the expression of
    all the annotated genes
    """

    #1. Load data and prepare a dictionary
    files = ['/home/smiravet/crg/mycorepo/minusTSSTTS.csv', '/home/smiravet/crg/mycorepo/plusTSSTTS.csv']
    expression = load_expression('./datasets/dsspilesmpn.txt')
    meanstd   = {}
    expressed = []

    #2. Assign mean and standard deviation to all the annotations
    for i in [0,1]:
        with open(files[i], 'r') as fi:
            for line in fi:
                genename, start, end = line.strip().split()
                start, end = int(start), int(end)
                if re.match('^[A-Z]{3}[0-9]{3}$', genename):
                    meanstd[genename] = [np.mean(expression.iloc[start:end,i]), np.std(expression.iloc[start:end,i])]
                    expressed        += list(expression.iloc[start:end, i])


    #3. Isolate N% most expressed genes, sort the dictionary and extract
    # the value mean-1SD for the lowest expressed in that group:

    if percentage > 0:
        number_genes = int(round((len(meanstd)*percentage)/100.0))
        top_genes = sorted(meanstd.items(), key=lambda i: i[1][0], reverse=True)[:number_genes]
        lowestinhighest = top_genes[-1]

        no_expression_1SD = lowestinhighest[1][0]-lowestinhighest[1][1]
        no_expression_2SD = lowestinhighest[1][0]-(lowestinhighest[1][1]*2)
    else:
        # Obtain the value for general expression
        tmean = np.mean(expressed)
        tstdv = np.std(expressed)

        no_expression_1SD = tmean
        no_expression_2SD = tmean - tstdv

    return no_expression_1SD, no_expression_2SD


#####################
#      CLASSES      #
#####################


#####################
#     EXECUTION     #
#####################

if __name__ == "__main__":
    print(no_expression_guesser(0))
    # merge_expression_experiments('/home/smiravet/Dropbox/mycorepo/systems_biology/Transcriptome/Rawdata/')

