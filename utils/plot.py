#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Plotting
Litt/Perry Group
Columbia University Mathematics REU 2018
'''

from functools import reduce
from matplotlib.collections import LineCollection

import matplotlib.pyplot as plt
import numpy as np
import csv
import colorsys
import fractions

# Functions #

# Utilities

def lcm(a, b):
    return abs(a * b) / fractions.gcd(a, b) if a and b else 0

def get_i_mod_p(y, p, i):
    return [n for n in y if (n - 1) % p == 0]

def plot_lines(list_of_vertices, savefile=None):
    list_of_lines = [[list_of_vertices[i], list_of_vertices[i+1]] for i in range(len(list_of_vertices) - 1)]
    #print list_of_vertices
    #print list_of_lines
    lc = LineCollection(list_of_lines)
    #print lc
    fig, ax = plt.subplots()
    xcors = [vertex[0] for vertex in list_of_vertices]
    ycors = [vertex[1] for vertex in list_of_vertices]
    ax.set_xlim(-1, max(xcors)+1)
    ax.set_ylim(min(-1, min(ycors)), max(ycors)+1)
    ax.add_collection(lc)
    plt.scatter(xcors, ycors)
    plt.grid(True)
    if savefile is not None:
        fig.savefig(savefile)
    else:
        plt.show()

def plot_multiple_lines(dict, savefile=None):
    fig, ax = plt.subplots()
    xmax = 1
    ymin = 0
    ymax = 1
    legend = [[], []]
    hue_step = 1.0/(len(dict) + 1)
    #print(hue_step)
    hue = 0
    for prime, vertices in dict.items():
        list_of_lines = [[vertices[i], vertices[i+1]] for i in range(len(vertices) - 1)]
        #print(list_of_lines)
        color = colorsys.hsv_to_rgb(hue, 1, 1)
        hue += hue_step
        lc = LineCollection(list_of_lines, color=color)
        legend[0].append(lc)
        legend[1].append(prime)
        xcors = [vertex[0] for vertex in vertices]
        xmax = max(xmax, max(xcors))
        ycors = [vertex[1] for vertex in vertices]
        ymin = min(ymin, min(ycors))
        ymax = max(ymax, max(ycors))
        ax.add_collection(lc)
    ax.set_xlim(0, xmax*1.3)
    ax.set_ylim(ymin, ymax)
    ax.legend(tuple(legend[0]), tuple(legend[1]))
    if savefile is not None:
        fig.savefig(savefile)
    else:
        plt.show()

def complex_roots(min_prime=2, max_prime=100, 
    min_exponent=3, max_exponent=10):
    pass

def prime_mod_gcd_exp(p, n):
    n, p = map(int, n.split('_')), np.array(int(p))
    list_gcd = np.array(reduce(gcd, n))
    warning = (p + 1) % list_gcd == 0

    return p[warning], list_gcd[warning], p[~warning], list_gcd[~warning]


# Main

def main(filepath='data/surface_search_over_primes.1528487580.supersingular.csv'):
    n, p = [], []

    with open(filepath,'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            p.append(int(row[0]))
            n.append(map(int, row[1].split("_")))

    list_gcd = np.array([reduce(gcd, i) for i in n])
    p = np.array(p)

    warning = (p + 1) % list_gcd == 0

    g = (sns.jointplot(p[warning], list_gcd[warning], 
        kind="scatter", stat_func=None,label='mod gcd', 
        color="r").set_axis_labels("Primes", "Exponents GCD"))
    h = (sns.jointplot(p[~warning], list_gcd[~warning], 
        kind="scatter", stat_func=None,label='mod gcd', 
        color="b").set_axis_labels("Primes", "Exponents GCD"))

    plt.show()

if __name__ == '__main__':
    prime_mod_gcd_exp(7, '2_2_4_7')
