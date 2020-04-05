#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Newton Polygons
Litt/Perry Group
Columbia University Mathematics REU 2018
'''

# Imports #

from sage.all import *
from rational_function import *
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
from point_counting import *
from rational_function import *
import itertools
from glob import glob
import os

def primesfrom2to(n):
    '''
    https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
    Input n>=6, Returns a array of primes, 2 <= p < n
    '''
    sieve = np.ones(n/3 + (n%6==2), dtype=np.bool)
    sieve[0] = False
    for i in xrange(int(n**0.5)/3+1):
        if sieve[i]:
            k=3*i+1|1
            sieve[((k*k)/3)::2*k] = False
            sieve[(k*k+4*k-2*k*(i&1))/3::2*k] = False
    return list(np.r_[2,3,((3*np.nonzero(sieve)[0]+1)|1)])


def get_intersect(a1, a2, b1, b2):
    """ 
    Returns the point of intersection of the lines passing through a2,a1 and b2,b1.
    a1: [x, y] a point on the first line
    a2: [x, y] another point on the first line
    b1: [x, y] a point on the second line
    b2: [x, y] another point on the second line
    """
    s = np.vstack([a1,a2,b1,b2])        # s for stacked
    h = np.hstack((s, np.ones((4, 1)))) # h for homogeneous
    l1 = np.cross(h[0], h[1])           # get first line
    l2 = np.cross(h[2], h[3])           # get second line
    x, y, z = np.cross(l1, l2)          # point of intersection
    if z == 0:                          # lines are parallel
        return (float('inf'), float('inf'))
    return (x/z, y/z)

# Functions #

def surface_newton_polygon(rational_function, p, q):
    poly = rational_function[1]
    R = QQ['T']
    t = R.gen()
    poly = R(poly/(1 - t)/(1 - q**2 * t))
    return poly.change_ring(Qp(p)).newton_polygon().vertices()

def default_newton_polygon(rational_function, p, q):
    poly = rational_function[1]
    return poly.change_ring(Qp(p)).newton_polygon().vertices()

def full_plot(max_prime=100, min_exponent=3, 
    max_exponent=10, m=60, affine=True):
    used_exponents = []
    coeffs = [Integer(1)]*4
    exponent_sets = map(list, sorted(itertools.product(*[range(min_exponent, max_exponent + 1)]*len(coeffs))))
    for p in primesfrom2to(max_prime):
        for exponents in exponent_sets:
            if os.path.isfile('data/report/newton_polygons/{}/{}_{}.png'.format('affine' if affine is True else 'weighted_projective', p, '_'.join(map(str, exponents)))):
                print 'Skipping Preexisting Case Prime {}, Exponents {}'.format(p, '_'.join(map(str, exponents)))
                continue
            if exponents in used_exponents:
                print 'Skipping Used Exponents Prime {}, Exponents {}'.format(p, '_'.join(map(str, exponents)))
                continue
            used_exponents.append(exponents)
            try:
                q = p
                sols, N_vals = affine_solutions(len(coeffs) - 1, p, q, m, coeffs, exponents)
                rational_function = find_rational(sols, 4, 5)
                poly = rational_function[1]
                deg = poly.degree()
                R = QQ['T']
                t = R.gen()
                div = (1 - t)/(1 - q**2 * t)
                poly = (poly/div).numerator()
                poly = R(poly)
                lcoeff = poly.coeffs()[-1]
                k = valuation(lcoeff, p)
                points = [(0, 0), (deg, 2*k), (deg, k)]
                list_of_vertices = rational_function[1].change_ring(Qp(p)).newton_polygon().vertices()
                list_of_lines = [[list_of_vertices[i], list_of_vertices[i+1]] for i in range(len(list_of_vertices) - 1)]
                lc = LineCollection(list_of_lines)
                fig, ax = plt.subplots()
                xcors = [vertex[0] for vertex in list_of_vertices]
                ycors = [vertex[1] for vertex in list_of_vertices]
                # ax.set_xlim(-1, max(xcors)+1)
                # ax.set_ylim(min(-1, min(ycors)), max(ycors)+1)
                ax.add_collection(lc)
                ax.scatter([x[0] for x in points], [x[1] for x in points], color='r')
                plt.grid(True)
                plt.title('Prime {}, Exponents {}'.format(p, exponents))
                fig.savefig('data/report/newton_polygons/{}/{}_{}.png'.format('affine' if affine is True else 'weighted_projective', p, '_'.join(map(str, exponents))))
            except:
                print 'ERROR Prime {}, Exponents {}'.format(p, '_'.join(map(str, exponents)))

# Execution #

if __name__ == '__main__':
    full_plot()
