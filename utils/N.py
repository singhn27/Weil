#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
N
Litt/Perry Group
Columbia University Mathematics REU 2018
Navtej Singh <singhnav@umich.edu>
'''

from sage.all import *

# Functions #

def iterate_N(alpha_num, chars_of_coeffs, gauss_sums, q):
    iterate = reduce(lambda x, y: x * y, 
        [chars_of_coeffs[i][alpha_num[i]] * gauss_sums[i][alpha_num[i]] 
        for i in xrange(len(alpha_num))], 1)
    return iterate/q
