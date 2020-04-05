#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Cayley Table
Litt/Perry Group
Columbia University Mathematics REU 2018
Navtej Singh <singhnav@umich.edu>
'''

# Imports #

from sage.all import *
from sage.matrix.operation_table import OperationTable
import itertools

# Classes #

class CayleyTable(object):
    def __init__(self, order):
        self.order = order
        self.multable = OperationTable(GF(self.order), operator.mul)
        self.addtable = OperationTable(GF(self.order), operator.add)

    def mult(self, a, b):
        return self.multable.table()[a][b]

    def add(self, a, b):
        return self.addtable.table()[a][b]

    def quartics(self):
        return [self.mult(self.mult(self.mult(x, x), x), x) for x in range(self.order)]

    def X_4_4(self, S):
        return reduce(lambda x, y: x + y, [int(self.add(self.add(self.add(w, x), y), z) == 0) for (w, 
            x, y, z) in [tup for tup in itertools.product(S, S, S, S)]], 0)

# Functions #

def main(order):
    ct = CayleyTable(order)
    n_solns = ct.X_4_4(ct.quartics())
    print n_solns

# Execution #

if __name__ == '__main__':
    main(9)
