#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Main
Litt/Perry Group
Columbia University Mathematics REU 2018
'''

# Imports #

from sage.all          import *
from newton_polygons   import *
from plot              import plot_lines
from point_counting    import *
from rational_function import *
from supersingular     import *
from hodge             import *
from smooth            import *

# Functions #

# Main

def run_projective(coeffs, powers, p, q, m = 60):
    betti = -1
    if (len(coeffs) == 4):
        L = lcm(powers)
        betti = sum(hodge(L, [L/powers[i] for i in range(4)])) + 3
    if (betti > 0 and my_WP_smooth(coeffs, powers, p, q, 1)):
        sols, N_vals = projective_solutions(len(coeffs) - 1, p, q, 2*betti + 10, coeffs, powers)
        ratl_fn = find_rational(sols, 1, 0, betti + 1)
    else:
        sols, N_vals = projective_solutions(len(coeffs) - 1, p, q, m, coeffs, powers)
        ratl_fn = find_rational(sols, 4, 5)

    return sols, N_vals, ratl_fn, is_supersingular(ratl_fn[0], ratl_fn[1], q)

def run_affine(coeffs, powers, p, q, m = 60):
    sols, N_vals = affine_solutions(len(coeffs) - 1, p, q, m, coeffs, powers)
    ratl_fn = find_rational(sols, 4, 5)

    return sols, N_vals, ratl_fn, is_supersingular(ratl_fn[0], ratl_fn[1], q)

def main():
    p = Integer(raw_input("Prime [3]: ") or 3)
    q = Integer(raw_input("Order of Finite Field [3]: ") or 3)
    f = log(q, p)

    while not f in ZZ:
        print("q is not a power of p")
        p = Integer(raw_input("Prime [3]: ") or 3)
        q = Integer(raw_input("Order of Finite Field [3]: ") or 3)
        f = log(q, p)

    r = Integer(raw_input("Number of Terms [4]: ") or 4)
    coeffs, powers = [], []

    for i in range(r):
        coeffs.append(Integer(raw_input("Coefficient [1]: ") or 1))
        powers.append(Integer(raw_input("Power [4]: ") or 4))

    r -= 1

    ambient_space = raw_input("Weighted Projective (WP) or Affine Solutions (A) [WP] : ")
    m = Integer(raw_input("Maximum Computation Degree [60]: ") or 60)
    wp = True
    if ambient_space == "WP" or ambient_space == "":
        sols, N_vals, ratl_fn, is_sup =  run_projective(coeffs, powers, p, q, m)
    else:
        wp = False
        sols, N_vals, ratl_fn, is_sup =  run_affine(coeffs, powers, p, q, m)
    print(sols)
    if ratl_fn == (0,0):
        print "No Recurrence Found"
        return sols, None, None, None
    elif wp and r == 3:
        newton_poly = surface_newton_polygon(ratl_fn, p, q)
    else:
        newton_poly = "The definition of newtown polygon is unknown"
    return sols, ratl_fn, is_sup, newton_poly

# Execution #

if __name__ == '__main__':
    sols, ratl_fn, is_sup, newton_poly = main()
    print ratl_fn
    print is_sup
    print newton_poly
