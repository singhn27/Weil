#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Gauss Sum Factorization
Litt/Perry Group
Columbia University Mathematics REU 2018
'''

# Imports #

from sage.all                   import *
from affine_patch_calculations  import *

# Functions #

def s(n, p, q):
    n = int(n)
    p = int(p)
    q = int(q)
    n = n % (q - 1)
    S = 0
    while n:
        S = S + n % p
        n = n / p
    return S

def calc_exceptional_poly(p, r, coeffs, powers):
    PolyRing = UCF['T']
    t = PolyRing.gen()
    m = lcm(powers)
    f = Mod(p, m).multiplicative_order()
    q = p**f
    alphas = find_exceptional_alpha(p, r, coeffs, powers)
    poly = PolyRing(1)

    for alpha_mult in alphas:
        R = (-1)**(r) / Integer(q)
        for i, alpha in enumerate(alpha_mult[1]):
            R *= alpha_mult[0]
            R *= gauss_sum(alpha, p, f)
        poly *= (t*R + 1)
    return poly.change_ring(QQ)

# give default coeffs to return multiplicitiy
def find_exceptional_alpha(p, r, coeffs, powers):
    m = lcm(powers)
    f = Mod(p, m).multiplicative_order()
    q = p**f

    alphas = get_alpha_tuples(r, p, f, coeffs, powers)
    exceptional_alpha = []

    for alpha in alphas:
        for mu in Integers(m).list_of_elements_of_multiplicative_group():
            S = 0
            for i in range(len(alpha[1])):
                j_num = alpha[1][i] * m
                S += s((q - 1) * mu * j_num / m, p, q)
            if S != (p-1) * f * (r + 1)/2:
                exceptional_alpha.append(alpha)
                break
        #print alpha_num, S, (p-1)*f * (r + 1)/2
    return exceptional_alpha

def count_exceptional_alpha(p, r, powers):
    m = lcm(powers)
    f = Mod(p, m).multiplicative_order()
    q = p**f

    alphas = get_sorted_alpha_tuples(r, p, f, [1] * len(powers), powers)
    P = 0
    n = 0

    for alpha in alphas:
        is_root_of_unity = True
        for mu in Integers(m).list_of_elements_of_multiplicative_group():
            S = 0
            for i in range(len(alpha[1])):
                j_num = alpha[1][i] * m
                S += s((q - 1) * mu * j_num / m, p, q)
            if S != (p-1) * f * (r + 1)/2:
                is_root_of_unity = False
                break
        if is_root_of_unity:
            P += alpha[0]
        else:
            n += alpha[0]
        #print alpha_num, S, (p-1)*f * (r + 1)/2
    return n, P

def exist_exceptional_alpha(p, r, powers):
    m = lcm(powers)
    f = Mod(p, m).multiplicative_order()
    q = p**f

    alphas = get_sorted_alpha_tuples(r, p, f, [1] * len(powers), powers)

    for alpha in alphas:
        S = 0
        for i in range(len(alpha[1])):
            j_num = alpha[1][i] * m
            S += s((q - 1) * j_num / m, p, q)
        if S != (p-1) * f * (r + 1)/2:
            return True
    return False



def compute_newton_polygon(r, p, powers):
    m = lcm(powers)
    f = Mod(p, m).multiplicative_order()
    q = p**f
    D = dict()
    alphas = get_sorted_alpha_tuples(r, p, f, [1] * len(powers), powers)
    for alpha in alphas:
        S = 0
        for a in alpha[1]:
            j_num = a * m
            S += s((q - 1) * j_num / m, p, q)
        if S not in D:
            D[S] = alpha[0]
        else:
            D[S] += alpha[0]

    point = (0,0)
    points = []
    points.append(copy(point))
    for i in sorted(D):
        slope = Integer(i) / Integer((p - 1) * f) - 1
        point = (point[0]  + D[i], point[1] + D[i] * slope)
        points.append(copy(point))

    return points
