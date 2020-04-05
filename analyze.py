#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Navtej Singh
MATH 575 UMich Ann Arbor
Adapted Freely from my other project 
https://www.github.com/singhn27/REU2018
'''

# Imports #

from sage.all  import *
from variety   import *
from tests     import primesfrom2to

# Functions #


def factorization_over_cyclotomic(V):
    n = V.degree
    p = V.characteristic
    z_fn = zeta_function(V)
    V.print_V()
    formated_numerator = print_rational(z_fn[0], V.order)
    formated_denominator = print_rational(z_fn[1], V.order)
    ratl_fn_string = "{}\n-------------------\n{}".format(formated_numerator, formated_denominator)
    print "Is this variety smooth? {}\n".format(V.is_smooth())
    print "Number of Solutions: {}\n".format(V.num_points())
    print "Zeta Function: \n {}\n".format(ratl_fn_string)
    poly = z_fn[1]
    if poly == 0:
        print "ERROR: zeta function not found"
        return None
    # in the case of a projective surface we ignore the top and bottom homology
    if isinstance(V, ProjectiveVariety) and V.r == 3:
        R = QQ['T']
        t = R.gen()
        poly = R(poly/(1 - t)/(1 - V.order**2 * t))
    newton_poly = poly.change_ring(Qp(V.characteristic)).newton_polygon().vertices()
    print "Newton Polygon: {}\n".format(newton_poly)
    print "Is this variety supersingular? {}\n".format(is_supersingular(z_fn[0], z_fn[1], V.order))
    print "Factor Denominator:"
    poly = z_fn[1].change_ring(CyclotomicField(n))
    roots = [(CC(r), m, p) for r,m in poly.roots()]

    print "Factorization Over Q(zeta{}): {}".format(n, poly.factor())
    poly = z_fn[1].change_ring(GF(p))
    if poly != 0:
        fact = poly.factor()
    else:
        fact = 0
    print "Factorization Over F_{}: {}\n".format(p, fact)
    print "Factor Cyclotomic Polynomial of degree {}:".format(n)
    poly =  cyclotomic_polynomial(n).change_ring(GF(p))
    if poly != 0:
        fact = poly.factor()
    else:
        fact = 0
    print "Factorization Over F_{}: {}".format(p, fact)
    return roots


def go():
    complex_points = []
    for p in primesfrom2to(167):
        V = ProjectiveVariety(Integer(p), DefiningEquation([1,1], [3,1]))
        complex_points += factorization_over_cyclotomic(V)
        print "\n\n\n"
    pts = [r for r,m,p in complex_points]
    return pts

if __name__ == '__main__':
    pts = go()
