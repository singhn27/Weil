#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Reductions
Litt/Perry Group
Columbia University Mathematics REU 2018
Navtej Singh <singhnav@umich.edu>
'''

# Imports #

from sage.all import *

# Functions #

# reduces powers to isomorphic variety of lower degree

def exp_reduce(powers):
    return [ gcd(powers[i], lcm([n_j for (j, n_j) in enumerate(powers) if j != i])) for i in range(len(powers)) ]
# resolves singularities

def exp_resolve(powers, p):
    return [ powers[i] / p**valuation(powers[i], p) for i in range(len(powers)) ]
# reduces exponents to a base variety which has equal zeta function and a rational map to the origional

def exp_base(powers, p):
    return exp_reduce(exp_resolve(powers, p))
