#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Point Counting
Litt/Perry Group
Columbia University Mathematics REU 2018
'''

# Imports #

from hodge             import *
from rational_function import *
from sage.all          import *
from supersingular     import *
from utils.alpha       import *
from utils.characters  import *
from utils.d           import *
from utils.N           import *
from utils.bin_lists   import *

# Functions #

def patch_sols(r, p, q, m, coeffs, powers, sums, chars, k):
    N_vals = []
    d = find_d(powers , q**k) # Compute the list of d_i
    N = Integer(0)
    alpha_num = [1 for i in d] # Store the numerator of alpha_i for each iteration
    counter = 0 # Number of alpha tuples tested
    end_counter = prod([i - 1 for i in d]) # Total number of alpha tuples
    while (counter < end_counter): # While some alpha tuples are untested
        if sum([alpha_num[i]/d[i] for i in range(0, len(d))]) in ZZ: # Check if the sum of alpha_i is an integer
            N += iterate_N(alpha_num, chars, sums, q**k)
        alpha_num = increment_alpha(alpha_num, d)
        counter += 1
    if type(N) != type(Integer(0)):
        N = Integer(N.to_cyclotomic_field())
        N += Integer(0)
    N_vals.append(N)
    N *= Integer(q**k - 1)
    N += Integer(q**(r*k))
    return N, N_vals

def find_X(u, p, q, m, coeffs, powers, sums, chars, k):
    if sum(u) == 0:
        return 1,0
    new_coeffs = [coeffs[i] for i in range(len(coeffs)) if u.is_one(i)]
    new_powers = [powers[i] for i in range(len(powers)) if u.is_one(i)]
    new_sums = [sums[i] for i in range(len(powers)) if u.is_one(i)]
    new_chars = [chars[i] for i in range(len(powers)) if u.is_one(i)]
    return patch_sols(len(new_coeffs) - 1, p, q, m, new_coeffs, new_powers, new_sums, new_chars, k)

def find_d_bin(u, powers):
    L = lcm(powers)
    return gcd([L/powers[i] for i in range(len(powers)) if u.is_one(i)])

def projective_solutions(r, p, q, m, coeffs, powers):
    r = Integer(r)
    p = Integer(p)
    q = Integer(q)
    m = Integer(m)
    coeffs = [Integer(coeff) for coeff in coeffs]
    powers = [Integer(power) for power in powers]
    N_sols, N_values = [], []

    f = log(q, p)
    for k in range(1, m + 1): # Calculate the number of points on X(F_{p^k}) for k up to m
        gauss_sums = get_sums(p, f, k, powers, r) # Compute 2-dim matrix of gauss sums
        chars_of_coeffs = get_chars(p, f, k, powers, coeffs, r) # Compute 2-dim matrix of chi_{alpha_i}(a_i^{-1})

        u = all_one(r)
        X = dict()

        for t in u.all_less_than():
            X_val, N_vals = find_X(t, p, q, m, coeffs, powers, gauss_sums, chars_of_coeffs, k)
            X[t.to_int()] = X_val
            N_values.extend(N_vals)
        C = dict()

        for t in u.all_less_than():
            C[t.to_int()] = alt_sum_over_bins(X, t) + (-1)**(t.sum())
        N = 0

        for t in u.all_less_than():
            N += C[t.to_int()] * Integer(gcd(find_d_bin(t, powers), q**k - 1)) / Integer(q**k - 1)

        N_sols.append(N)

    return N_sols, N_values

def affine_solutions(r, p, q, m, coeffs, powers):
    r = Integer(r)
    p = Integer(p)
    q = Integer(q)
    m = Integer(m)
    coeffs = [Integer(coeff) for coeff in coeffs]
    powers = [Integer(power) for power in powers]
    N_sols, N_values = [], []

    f = log(q, p)

    for k in range(1, m + 1): # Calculate the number of points on X(F_{p^k}) for k up to m
        gauss_sums = get_sums(p, f, k, powers, r) # Compute 2-dim matrix of gauss sums
        chars_of_coeffs = get_chars(p, f, k, powers, coeffs, r) # Compute 2-dim matrix of chi_{alpha_i}(a_i^{-1})
        N, N_vals = patch_sols(r, p, q, m, coeffs, powers, gauss_sums, chars_of_coeffs, k)
        N_sols.append(N)
        N_values.extend(N_vals)

    return N_sols, N_values
