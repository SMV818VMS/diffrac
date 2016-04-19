# Diffraction

## Goal

Diffraction includes several function to work and perform differential analysis over pile ups coming from non fractionated RNASeq analysis.


## finddrops

This script looks for drops in the expression using the following algorithm:

1. _finddrops_ accepts a two columns files ([position, expression])
2. It takes the expression profile and an annotation file to compute the average distance between annotations.
3. It works with sliding windows composed by a expression window (100 by default) and a no\_exp\_window (average distance between annotations/4)
4. The drop will only be considered if:
  - Window of expression does not include any 0.0 and the last nucleotide is greater than 0.0
  - Window of no expression has an average smaller than a threshold of expression (0.0 by default)
5. The outcome is a dictionary including: \{identifier:\[window start, window end, last position with expression, stdsc, maxsc, dropsc\]\}
  - As you can notice we see three different scores: stdsc, maxsc and dropsc:
    - stdsc is the standard deviation of the expression window
    - maxsc is the maximum difference between contiguous positions
    - dropsc is the first derivative between the last position and the first nucleotide in the no expression window

The interpretation for this values is:

  - high stdsc, low maxsc and maxsc != dropsc ==> decay position: The change in expression decays gradually but with no big changes, this hardens the match between maximum change (minimum derivative) and the drop in the last position with expression
  - low stdsc, high maxsc and maxsc == dropsc ==> sharp position: The change in expression is abrupt and with a big difference, this makes easier the match between minimum derivative and the index position.

The goal is to have a simple rule to discriminate between:


**Decay drop**
![decay drop](figures/decay)

**Sharp drop**
![sharp drop](figures/sharp)


