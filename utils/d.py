#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
d
Litt/Perry Group
Columbia University Mathematics REU 2018
Ming Jing <mjing@fordham.edu>
'''

# Imports #

from sage.all import *

# Functions #

def find_d(n, q):
    return [Integer(gcd(n_i, q-1)) for n_i in n]
