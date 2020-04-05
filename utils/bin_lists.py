#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Weighted Projective
Litt/Perry Group
Columbia University Mathematics REU 2018
Benjamin Church <bvc2105@columbia.edu>
'''

# Classes #

class bin_list:

    def __init__(self):
        self.data = []
        self.length = 0

    def __init__(self, __data__):
        self.data = list(__data__)
        self.length = len(self.data)

    def to_int(self):
        s = 0
        for i in range(self.length):
            if self.data[i] == 1:
                s += 2**i
        return s

    def __hash__(self):
        return self.to_int()

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def bin_print(self):
        for i in self.data:
            print i,
        print "Int: ", self.to_int()

    def is_zero(self):
        for i in self.data:
            if i != 0:
                return False
        return True
        
    def is_less_than(self, r):
        if self.length != r.length:
            return False
        else:
            for i in range(r.length):
                if r.data[i] == 0 and self.data[i] != 0:
                    return False
            return True

    def all_less_than(self):
        if self.is_zero():
            return {}
        else:
            ls = {self}
            for i in range(self.length):
                if self.data[i] != 0:
                    new_dat = list(self.data)
                    new_dat[i] = 0
                    new_bin = bin_list(new_dat)
                    ls = ls.union(new_bin.all_less_than())
            return ls

    def is_one(self, i):
        return (self.data[i] == 1)

    def sum(self):
        return sum(self.data)

def all_one(r):
    return bin_list([1]*(r+1))

def print_bins(bs):
    for b in bs:
        b.bin_print()


def alt_sum_over_bins(X, u):
    s = 0
    for t in u.all_less_than():
        s += (-1)**(u.sum() - t.sum()) * X[t.to_int()]
    return s

"""
u = bin_list([1,1,1,1])
u.bin_print()
print ""
print_bins(u.all_less_than())
print ""

C = dict()

for t in u.all_less_than():
    C[t.to_int()] = 1

print C

print sum_over_bins(C, bin_list([1,1,1,1])), sum_over(C, bin_list([0,1,1,1])), sum_over(C, bin_list([0,0,1,1]))

print u.__hash__()
"""
