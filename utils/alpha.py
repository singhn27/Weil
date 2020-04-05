#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Alpha Calculations
Litt/Perry Group
Columbia University Mathematics REU 2018
'''

# Imports #

from sage.all   import *
from d          import *
from bin_lists  import *
from characters import *

# Functions #

def find_d_bin(u, powers):
    L = lcm(powers)
    return gcd([L/powers[i] for i in range(len(powers)) if u.is_one(i)])


def maximal_alpha_num(r, powers):
    u = all_one(r)
    S = 0
    for t in u.all_less_than():
        t_prod = 1
        t_powers = list()
        for i in range(r+1):
            if t.is_one(i):
                t_prod *= powers[i]
                t_powers.append(powers[i])
        S += (-1)**(r+1 - t.sum()) * t_prod / lcm(t_powers)
    # include the all zero case
    S += (-1)**(r+1)
    return S

# gets the next valid alpha

def is_alpha_integer(alpha, d):
    return sum([Integer(alpha[i]) / Integer(d[i]) for i in range(0, len(d))]) in ZZ

def make_alpha_integer(alpha, d):
    r = len(d) - 1
    s = sum([Integer(alpha[i]) / Integer(d[i]) for i in range(0, r)])
    needed = (1 - frac(s))*d[r]
    # if the leftover bit is greater than zero and of the form a/d[r]
    if needed > 0 and needed < d[r] and needed in ZZ:
        alpha[r] = needed
        return alpha
    else:
        return None


def next_alpha(alpha, d):
    r = len(d) - 1
    while [alpha[i] for i in range(0,r)] != [d[i] - 1 for i in range(0, r)]:
        i = 0
        while i < r:
            if alpha[i] < d[i] - 1:
                alpha[i] += 1
                break
            alpha[i] = 1
            i += 1
        new_alpha = make_alpha_integer(alpha, d)
        if new_alpha is not None:
            return new_alpha
    return None

# Legacy function increments alpha by one

def increment_alpha(alpha, d):
    i = 0
    while i < len(d):
        if alpha[i] < d[i] - 1:
            alpha[i] += 1
            return alpha
        alpha[i] = 1
        i += 1
    return alpha

def tuple_cmp(lista, listb):
    for i in range(len(lista[1])):
        if lista[1][i] > listb[1][i]:
            return 1
        if lista[1][i] < listb[1][i]:
            return -1
    return 0

def equals_upto(lista, listb, index):
    return all([lista[i] == listb[i] for i in range(index)])


def get_unsorted_alphas(r, p, k, powers):
    d = find_d(powers , p**k) # Compute the list of d_i

    alpha_num = [1 for i in d] # Store the numerator of alpha_i for each iteration
    new_alpha = make_alpha_integer(alpha_num, d)
    if new_alpha is not None:
        alpha_num = new_alpha

    counter = 0 # Number of alpha tuples tested
    end_counter = maximal_alpha_num(r, d) # Total number of alpha tuples

    #Get possible alphas
    alphas = []
    while (counter < end_counter): # While some alpha tuples are untested
        if is_alpha_integer(alpha_num, d): # Check if the sum of alpha_i is an integer
            #Add the (product of chars, list of alphas) tuple to alpha_tuples
            alphas.append([Integer(alpha_num[i])/Integer(d[i]) for i in range(r + 1)])
        alpha_num = next_alpha(alpha_num, d)
        if alpha_num == None:
            break
        counter += 1
    return alphas

def get_alpha_tuples(r, p, k, coeffs, powers):
    alphas = get_unsorted_alphas(r, p, k, powers) # Get all possible alpha lists
    # Calculate characters for each
    alpha_tuples = []
    for alpha in alphas:
        char_prod = prod([char(alpha[i], coeffs[i], p, k) for i in range(0, len(coeffs))])
        alpha_tuples.append((conjugate(char_prod), alpha))
    return alpha_tuples

# give default coeffs to return multiplicitiy
def get_sorted_alpha_tuples(r, p, k, coeffs, powers):
    alpha_tuples = get_alpha_tuples(r, p, k, coeffs, powers) # get the unsorted list of alpha tuples
    alpha_tuples = [(alpha[0], sorted(alpha[1])) for alpha in alpha_tuples] #Sort the list of alphas in each tuple
    alpha_tuples = sorted(alpha_tuples, cmp=tuple_cmp) #Sort the tuples by the alpha component

    #Combine tuples with equal alpha component by summing their 0-place
    alpha_mults = []
    i = 0
    while i < len(alpha_tuples):
        j = i
        char_sum = Integer(0)
        while j < len(alpha_tuples) and alpha_tuples[i][1] == alpha_tuples[j][1]:
            char_sum += alpha_tuples[j][0]
            j += 1
        alpha_mults.append((char_sum, alpha_tuples[i][1]))
        i = j
    return alpha_mults

# Legacy Functions

def return_alphas(r, p, k, coeffs, powers):
    d = find_d(powers , p**k) # Compute the list of d_i
    alpha_num = [1 for i in d] # Store the numerator of alpha_i for each iteration
    counter = 0 # Number of alpha tuples tested
    end_counter = prod([i - 1 for i in d]) # Total number of alpha tuples

    #Get possible alpha tuples
    alpha_tuples = []
    while (counter < end_counter): # While some alpha tuples are untested
        if sum([Integer(alpha_num[i])/Integer(d[i]) for i in range(0, len(d))]) in ZZ: # Check if the sum of alpha_i is an integer
            #Add the (product of chars, list of alphas) tuple to alpha_tuples
            alpha_tuples.append((char(Integer(alpha_num[i])/Integer(d[i]), coeffs[i], p, k), [Integer(alpha_num[i])/Integer(d[i]) for i in range(r + 1)]))
        alpha_num = increment_alpha(alpha_num, d)
        counter += 1
    alpha_tuples = [(alpha[0], sorted(alpha[1])) for alpha in alpha_tuples] #Sort the list of alphas in each tuple
    alpha_tuples = sorted(alpha_tuples, cmp=tuple_cmp) #Sort the tuples by the alpha component

    #Combine tuples with equal alpha component by summing their 0-place
    alpha_mults = []
    i = 0
    while i < len(alpha_tuples):
        j = i
        char_sum = Integer(0)
        while j < len(alpha_tuples) and alpha_tuples[i][1] == alpha_tuples[j][1]:
            char_sum += alpha_tuples[j][0]
            j += 1
        alpha_mults.append((char_sum, alpha_tuples[i][1]))
        i = j
    return alpha_mults

def test_alphas(r, p, k, coeffs, powers):
    d = find_d(powers , p**k) # Compute the list of d_i
    alpha_num = [1 for i in d] # Store the numerator of alpha_i for each iteration
    counter = 0 # Number of alpha tuples tested
    end_counter = maximal_alpha_num(r, d) # Total number of alpha tuples

    #Get possible alpha tuples
    alpha_tuples = []
    while (counter < end_counter): # While some alpha tuples are untested
        if sum([Integer(alpha_num[i])/Integer(d[i]) for i in range(0, len(d))]) in ZZ: # Check if the sum of alpha_i is an integer
            #Add the (product of chars, list of alphas) tuple to alpha_tuples
            char_prod = prod([char(Integer(alpha_num[i])/Integer(d[i]), coeffs[i], p, k) for i in range(0, len(d))])
            alpha_tuples.append((conjugate(char_prod), [Integer(alpha_num[i])/Integer(d[i]) for i in range(r + 1)]))
        alpha_num = next_alpha(alpha_num, d)
        counter += 1
    alpha_tuples = [(alpha[0], sorted(alpha[1])) for alpha in alpha_tuples] #Sort the list of alphas in each tuple
    alpha_tuples = sorted(alpha_tuples, cmp=tuple_cmp) #Sort the tuples by the alpha component

    #Combine tuples with equal alpha component by summing their 0-place
    alpha_mults = []
    i = 0
    while i < len(alpha_tuples):
        j = i
        char_sum = Integer(0)
        while j < len(alpha_tuples) and alpha_tuples[i][1] == alpha_tuples[j][1]:
            char_sum += alpha_tuples[j][0]
            j += 1
        alpha_mults.append((char_sum, alpha_tuples[i][1]))
        i = j
    return alpha_mults
