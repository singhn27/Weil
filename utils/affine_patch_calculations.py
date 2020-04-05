#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Point Counting
Litt/Perry Group
Columbia University Mathematics REU 2018
'''

# Imports #

from sage.all             import *
from alpha                import *
from characters           import *
from d                    import *
from bin_lists            import *

# Functions #

def compute_N_val(alpha_mult, r, p, k, start, index):
    N = Integer(0)
    if index == r + 1:
        return alpha_mult[start][0]
    j = 0
    while start + j < len(alpha_mult) and equals_upto(alpha_mult[start][1], alpha_mult[start + j][1], index):
        i = 0
        N += gauss_sum(alpha_mult[start + j][1][index], p, k) * compute_N_val(alpha_mult, r, p, k, start + j, index + 1)
        while start + i + j < len(alpha_mult) and equals_upto(alpha_mult[start + j][1], alpha_mult[start + j + i][1], index + 1):
            i += 1
        j += i
    return N



def patch_sols(r, p, k, coeffs, powers):
    alpha_mults = return_alphas(r, p, k, coeffs, powers)
    #Compute N
    N = Integer(0)
    N += compute_N_val(alpha_mults, r, p, k, 0, 0)
    N += Integer(0)
    if type(N) != type(Integer(0)): #ensure types are correct
        N = Integer(N.to_cyclotomic_field())
        N += Integer(0)
    N /= Integer(p**k)
    N *= Integer(p**k - 1)
    N += Integer(p**(r*k))
    return N


'''
from time import time
t = time()
print len(test_alphas(3, 97, 1, [1,1,1,1], [96,96,96,96]))
print time() - t
t = time()
print len(return_alphas(3, 97, 1, [1,1,1,1], [96,96,96,96]))
print time() - t
'''
#patch_sols(3, 7, 1, [1,1,1,1], [7,7,7,7])
# print get_alpha_nums_with_multiplicity(3, 5, 2, [4,4,4,10])
# print get_alpha_nums_with_multiplicity(3, 5, 2, [4,4,4,2])

#print test_alphas(5, [4,4,4,16,56,1244]), maximal_alpha_num(5, [4,4,4,16,56,1244])
