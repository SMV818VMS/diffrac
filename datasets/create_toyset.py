#!/usr/bin/env python

import random


# The idea is to create a toy dataset including the two types of termination
# we can see in the RNA seq profiles to test two different things:

# 1. The size of the windows required to detect the signals
# 2. The whole approach


# We will put the information in the toyset file
fo = open('./toyset.txt','w')

# Define the properties of the population and the position for 2 types of
# termination, decay and sharp. We define the list of positions to be able to test they work:

sharp = [2000, 4000, 6000, 8000]
decay = [1000, 3000, 5000, 7000]
window = 100

expression =  



for i in list(range(1, 100)):

