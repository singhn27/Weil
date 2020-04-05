#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Main
Litt/Perry Group
Columbia University Mathematics REU 2018
Navtej Singh <singhnav@umich.edu>
'''

# Imports #

import cProfile
import gc
import itertools
import signal

import numpy as np

from glob              import glob
from main              import run_affine, run_projective
from newton_polygons   import *
from plot              import plot_lines, plot_multiple_lines
from rational_function import *
from sage.all          import *
from smooth            import my_WP_smooth
from supersingular     import *
from time              import time
from variety           import AffineVariety, DefiningEquation, ProjectiveVariety


# Functions #

# Timeout

class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException

signal.signal(signal.SIGALRM, timeout_handler)

# Math

def primesfrom2to(n):
    '''
    https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
    Input n>=6, Returns a array of primes, 2 <= p < n
    '''
    sieve = np.ones(n/3 + (n%6==2), dtype=np.bool)
    sieve[0] = False
    for i in xrange(int(n**0.5) / 3 + 1):
        if sieve[i]:
            k=3*i+1|1
            sieve[((k*k)/3)::2*k] = False
            sieve[(k*k+4*k-2*k*(i&1))/3::2*k] = False
    return [Integer(prime) for prime in list(np.r_[2,3,((3*np.nonzero(sieve)[0]+1)|1)])]

def primitive_root(a, b):
    return Mod(a, b).multiplicative_order() == b - 1

# Tests

def multiple_newton_polygon_plot(min_prime=2, max_prime=2,
    min_exponent=3, max_exponent=4, affine=True):
    coeffs, np_vertices = [Integer(1)] * 4, []
    exponent_sets = map(list, sorted(itertools.product(*[range(min_exponent, max_exponent + 1)]*len(coeffs))))
    primes = map(lambda p: Integer(p), filter(lambda p: p >= min_prime, primesfrom2to(int(max_prime)+1)))
    for exponent_set in exponent_sets:
        for prime in primes:
            solutions, N_values, ratl_fn, supersingularity = run_affine(coeffs, exponent_set, prime, prime) if affine else run_projective(coeffs, exponent_set, prime, prime)
            np_vertices.append(default_newton_polygon(ratl_fn, prime, prime))
    [plot_lines(vertex_set) for vertex_set in np_vertices]

# Conditions

def condition_gcd_lcm_exp(exponent_set):
    return all([gcd(n_i, lcm([n_j for (j, n_j) in enumerate(exponent_set) if j != i])) == n_i for (i,
        n_i) in enumerate(exponent_set)])

def condition_increasing_exp(exponent_set):
    for i in range(0, len(exponent_set) - 1):
        if exponent_set[i] > exponent_set[i+1]:
            return False
    return True

def condition_sum_weights(exponent_set):
    L = lcm(exponent_set)
    return sum([L / n_i for n_i in exponent_set]) - L < 0

def condition_prime_congruent_lcm(prime, exponent_set):
    lcm_n_i, power = lcm(exponent_set), 1
    while prime ** power % lcm_n_i != 1:
        if power > lcm_n_i + 5:
            return None
        if prime ** power % lcm_n_i - lcm_n_i == -1:
            return True
        power += 1
    return False

def condition_coprime_pairs(exponent_set):
    for n in exponent_set:
        if len([ k for k in exponent_set if gcd(k,n) > 2 ]) > 2:
            return False
    return True


# Plotting

def variety_newton_polygon_plot(min_prime=2, max_prime=2, coeffs=[1,1,1,1], powers=[4,4,4,4], savefile=None):
    primes = map(lambda p: Integer(p), filter(lambda p: p >= min_prime, primesfrom2to(int(max_prime) + 1)))
    d = {}
    for prime in primes:
        solutions, N_values, ratl_fn, supersingularity = run_projective(coeffs, powers, prime, prime)
        d[prime] = surface_newton_polygon(ratl_fn, prime, prime)
    plot_multiple_lines(d, savefile)

# Data Generation

def faster_search_over_primes(affine=True, min_prime=2, max_prime=100, min_exponent=3, max_exponent=25, timeout=30, time_precision=4, profile=True):
    print '\nBeginning Run Case {} for Primes {} to {}, Exponents {} to {}, Timeout {}, Time Precision {}, {} Profiling\n'.format('Affine' if affine else 'WP',
        min_prime, max_prime, min_exponent, max_exponent, timeout, time_precision, 'With' if profile else 'Without')
    affine, cases, coeffs, errors, profile, start, used_exponent_sets = bool(affine), 0, [Integer(1)] * 4, 0, bool(profile), int(time()), []
    print 'Set Variables\n'
    exponent_sets = map(list, sorted(itertools.product(*[range(min_exponent, max_exponent + 1)] * len(coeffs))))
    print 'Generated Exponents\n'
    primes = map(lambda p: Integer(p), filter(lambda p: p >= min_prime, primesfrom2to(int(max_prime) + 1)))
    print 'Generated Primes\n'
    if profile is True:
        profiler = cProfile.Profile()
        print 'Instantiated Profiler\n'
    with open('data/{}/faster_search_over_primes.{}.total.csv'.format('affine' if affine else 'weighted_projective', start), 'a') as outfile:
        outfile.write('prime,exponents,smoothness,supersingularity,error,raw_speed,calls,callcount,reccallcount,inlinetime,totaltime\n')
        for exponent_set in exponent_sets:
            if exponent_set in used_exponent_sets:
                continue
            used_exponent_sets.append(exponent_set)
            for prime in primes:
                print 'Trying Prime {}, Exponents {}\n'.format(prime, '_'.join(map(str, exponent_set)))
                try:
                    tick = time()
                    variety = AffineVariety(prime, DefiningEquation(coeffs, exponent_set)) if affine is True else ProjectiveVariety(prime, DefiningEquation(coeffs, exponent_set))
                    smoothness = variety.is_smooth()
                    supersingularity = variety.is_supersingular()
                    raw_speed = time() - tick
                    if profile is True:
                        try:
                            stats = profiler.getstats()[-1]
                            line = '{},{},{},{},None,{},{},{},{},{},{}\n'.format(*map(str, [prime, '_'.join(map(str, exponent_set)), smoothness, supersingularity, round(raw_speed, time_precision),
                                stats.calls, stats.callcount, stats.reccallcount, round(stats.inlinetime, time_precision), round(stats.totaltime, time_precision)]))
                        except IndexError:
                            line = '{},{},{},{},None,{},N/A,N/A,N/A,N/A,N/A\n'.format(*map(str, [prime, '_'.join(map(str, exponent_set)), smoothness, supersingularity, round(raw_speed, time_precision)]))
                    else:
                        line = '{},{},{},{},None,{},N/A,N/A,N/A,N/A,N/A\n'.format(*map(str, [prime, '_'.join(map(str, exponent_set)), smoothness, supersingularity, round(raw_speed, time_precision)]))
                    print line
                    outfile.write(line)
                    cases += 1
                except TimeoutException as e:
                    if profile is True:
                        try:
                            stats = profiler.getstats()[-1]
                            line = '{},{},N/A,N/A,TimeoutException,N/A,{},{},{},{},{}\n'.format(*map(str, [prime, '_'.join(map(str, exponent_set)), stats.calls, stats.callcount, stats.reccallcount, stats.inlinetime, stats.totaltime]))
                        except IndexError:
                            line = '{},{},N/A,N/A,TimeoutException,N/A,N/A,N/A,N/A,N/A,N/A\n'.format(*map(str, [prime, '_'.join(map(str, exponent_set))]))
                    else:
                        line = '{},{},N/A,N/A,TimeoutException,N/A,N/A,N/A,N/A,N/A,N/A\n'.format(*map(str, [prime, '_'.join(map(str, exponent_set))]))
                    print line
                    outfile.write(line)
                    print '\nERROR (continuing...) TimeoutException\n'
                    errors += 1
                except Exception as e:
                    if profile is True:
                        try:
                            stats = profiler.getstats()[-1]
                            line = '{},{},N/A,N/A,{},N/A,{},{},{},{},{}\n'.format(*map(str, [prime, '_'.join(map(str, exponent_set)), e.message, stats.calls, stats.callcount, stats.reccallcount, stats.inlinetime, stats.totaltime]))
                        except IndexError:
                            line = '{},{},N/A,N/A,{},N/A,N/A,N/A,N/A,N/A,N/A\n'.format(*map(str, [prime, '_'.join(map(str, exponent_set)), e.message]))
                    else:
                        line = '{},{},N/A,N/A,{},N/A,N/A,N/A,N/A,N/A,N/A\n'.format(*map(str, [prime, '_'.join(map(str, exponent_set)), e.message]))
                    print line
                    outfile.write(line)
                    print '\nERROR (continuing...) {}\n'.format(e.message)
                    errors += 1
                    signal.alarm(0)
                finally:
                    gc.collect()
                    signal.alarm(0)
    end = int(time())
    print '\n\nFinished {} cases in {} seconds with {} errors\n\n'.format(cases, end - start, errors)

def loop_fixed_variety_over_primes(exponents, min_prime=2, max_prime=150, coeffs=[1, 1, 1, 1], affine=True):
    for prime in [p for p in primesfrom2to(max_prime) if p >= min_prime]:
        variety = AffineVariety(prime, DefiningEquation(coeffs, exponents)) if affine is True else ProjectiveVariety(prime, DefiningEquation(coeffs, exponents))
        print '{} Variety, Coefficients {}, Exponents {}, Supersingularity {}'.format('Affine' if affine is True else 'Projective', ' '.join(map(str, coeffs)), ' '.join(map(str, exponents)), variety.is_supersingular())

def my_test(min_prime=2, max_prime=150, min_exponent=20,
    max_exponent=70, timeout = 30, affine=True):
    print '\nBeginning Run Case {} for Primes {} to {}, Exponents {} to {}'.format('Affine' if affine else 'WP',
        min_prime, max_prime, min_exponent, max_exponent)
    cases, coeffs, errors, original_timeout, start, used_exponent_sets = [], \
        [Integer(1)] * 4, [], timeout, int(time()), []
    parameters = {}
    exponent_sets = map(list, sorted(itertools.product(*[range(min_exponent, max_exponent + 1)]*len(coeffs))))
    primes = map(lambda p: Integer(p), filter(lambda p: p >= min_prime, primesfrom2to(int(max_prime)+1)))
    for exponent_set in exponent_sets:
        exponent_set = '_'.join(map(str, exponent_set))
        for prime in primes:
            if exponent_set not in parameters:
                parameters[exponent_set] = [prime]
            else:
                parameters[exponent_set].append(prime)
    for exponent_set in exponent_sets:
        if not condition_gcd_lcm_exp(exponent_set) or not condition_increasing_exp(exponent_set) or condition_coprime_pairs(exponent_set):
            continue
        else:
            print "Checking {}".format(exponent_set)
            n = lcm(exponent_set)
            for p in primes:
                if n % p == 0:
                    continue
                else:
                    f = Mod(p, n).multiplicative_order()

                def_equation = DefiningEquation([1,1,1,1], exponent_set)
                V = AffineVariety(p, def_equation) if affine else ProjectiveVariety(p, def_equation)

                is_sm = V.is_smooth()
                is_SS = V.is_supersingular()

                if f > 0 and f % 2 == 0:
                    to_check = Mod(p**(f/2), n)
                    if to_check == Mod(-1, p) and not is_SS:
                        print "VIOLATION OF SHIODA\np = {}, powers = {}, f = {} and p**(f/2) ~ {} but Variety is NOT Supersingular\nIs Variety Smooth: {}".format(p, exponent_set, f, to_check, is_sm)
                    if to_check != -1 and is_SS:
                        print "VIOLATION OF SHIODA CONVERSE\np = {}, powers = {}, f = {} and p**(f/2) ~ {} but Variety is Supersingular\nIs Variety Smooth: {}".format(p, exponent_set, f, to_check, is_sm)
                elif is_SS:
                    print "VIOLATION OF SHIODA CONVERSE\np = {}, powers = {}, f = {} but Variety is Supersingular\nIs Variety Smooth: {}".format(p, exponent_set, f, is_sm)


def search_over_primes(min_prime=2, max_prime=100, min_exponent=3,
    max_exponent=10, timeout=30, parameters=None, affine=True):
    print '\nBeginning Run Case {} for Primes {} to {}, Exponents {} to {}, and Timeout {}\n'.format('Affine' if affine else 'WP',
        min_prime, max_prime, min_exponent, max_exponent, timeout)
    cases, coeffs, errors, original_timeout, start, used_exponent_sets = [], \
        [Integer(1)] * 4, [], timeout, int(time()), []
    if parameters is None:
        parameters = {}
        exponent_sets = map(list, sorted(itertools.product(*[range(min_exponent, max_exponent + 1)]*len(coeffs))))
        primes = map(lambda p: Integer(p), filter(lambda p: p >= min_prime, primesfrom2to(int(max_prime)+1)))
        for exponent_set in exponent_sets:
            exponent_set = '_'.join(map(str, exponent_set))
            for prime in primes:
                if exponent_set not in parameters:
                    parameters[exponent_set] = [prime]
                else:
                    parameters[exponent_set].append(prime)
    with open('data/{}/surface_search_over_primes.{}.total.csv'.format('affine' if affine else 'weighted_projective', start), 'a') as total_outfile:
        total_outfile.write('prime,exponents,numerator,denominator,supersingularity,smoothness,error,speed\n')
        for exponent_set, primes in sorted(parameters.iteritems(), key = lambda p: map(int, p[0].split('_'))):
            exponent_set = sorted(map(Integer, exponent_set.split('_')))
            if exponent_set in used_exponent_sets:
                continue
            used_exponent_sets.append(exponent_set)
            if not condition_gcd_lcm_exp(exponent_set):
                continue
            if not condition_sum_weights(exponent_set):
                timeout = 10
            for prime in primes:
                if condition_prime_congruent_lcm(prime, exponent_set):
                    timeout = 10
                signal.alarm(timeout)
                timeout = original_timeout
                try:
                    tick = time()
                    solutions, N_values, ratl_fn, supersingularity = run_affine(coeffs, exponent_set, prime, prime) if affine else run_projective(coeffs, exponent_set, prime, prime)
                    speed = time() - tick
                    if ratl_fn[0] == 0 and ratl_fn[1] == 0:
                        print 'Rational Function 0/0, Retrying with m = 120'
                        tick = time()
                        solutions, N_values, ratl_fn, supersingularity = run_affine(coeffs, exponent_set, prime, prime, 120) if affine else run_projective(coeffs, exponent_set, prime, prime, 120)
                        speed = time() - tick
                        if ratl_fn[0] == 0 and ratl_fn[1] == 0:
                            no_recurrence_error_line = '{},{},{}\n'.format(prime, '_'.join(map(str, exponent_set)), 'Rational Function 0/0')
                            print no_recurrence_error_line
                            continue
                    numerator = print_rational(ratl_fn[0], prime)
                    denominator = print_rational(ratl_fn[1], prime)
                    smoothness = my_WP_smooth(coeffs, exponent_set, prime, prime, 1)
                    line = '{},{},{},{},{}\n'.format(prime, '_'.join(map(str, exponent_set)), numerator, denominator, round(speed, 2))
                    full_line = '{},{},{},{},{},{},None,{}\n'.format(prime, '_'.join(map(str, exponent_set)), numerator, denominator, str(supersingularity), str(smoothness), round(speed, 2))
                    total_outfile.write(full_line)
                    if not supersingularity:
                        plot_lines(surface_newton_polygon(ratl_fn, prime, prime), savefile='data/{}/plots/{}.png'.format('affine' if affine else 'weighted_projective',
                            '_'.join(map(str, [prime, '_'.join(map(str, exponent_set))]))))
                    cases.append((prime, exponent_set))
                except TimeoutException as e:
                    full_line = '{},{},N/A,N/A,N/A,TimeoutException,N/A\n'.format(prime, '_'.join(map(str, exponent_set)))
                    total_outfile.write(full_line)
                    print '\nERROR (continuing...) TimeoutException\n'
                    errors.append(e)
                except Exception as e:
                    full_line = '{},{},N/A,N/A,N/A,{},N/A\n'.format(prime, '_'.join(map(str, exponent_set)), e.message)
                    total_outfile.write(full_line)
                    print '\nERROR (continuing...) {}\n'.format(e.message)
                    errors.append(e)
                    signal.alarm(0)
                finally:
                    gc.collect()
                    signal.alarm(0)
    end = int(time())
    print '\n\nFinished {} cases in {} seconds with {} errors\n\n'.format(len(cases), end - start, len(errors))

def rerun_search_over_primes(filepath='data/affine', wildcard='errors'):
    error_files = glob('{}/*{}*'.format(filepath, wildcard))
    parameters = {}
    for error_file in error_files:
        for index, line in enumerate(open(error_file, 'r')):
            try:
                prime, exponent_set, error_message = line.split(',')
                if exponent_set not in parameters:
                    parameters[exponent_set] = [prime]
                else:
                    parameters[exponent_set].append(prime)
            except ValueError:
                print 'Skipping Line {} in File {}'.format(index + 1, error_file)
    search_over_primes(parameters = parameters, affine = 'affine' in filepath)

def execute():
    affine = bool(int(raw_input('Affine (1) or Weighted Projective (0) [1]: ').strip().lower() or 1))
    min_prime = Integer(raw_input('Minimum Prime [2]: ') or 2)
    max_prime = Integer(raw_input('Maximum Prime [100]: ') or 100)
    min_exponent = Integer(raw_input('Minimum Exponent [3]: ') or 3)
    max_exponent = Integer(raw_input('Maximum Exponent [25]: ') or 25)
    timeout = int(raw_input('Timeout in Seconds [30]: ') or 30)
    time_precision = int(raw_input('Precision of Profile Reporting [6]: ') or 6)
    profile = int(raw_input('Profile (1) or Not (0) [0]: ') or 0)

    faster_search_over_primes(affine=affine, min_prime=min_prime, max_prime=max_prime, min_exponent=min_exponent,
        max_exponent=max_exponent, timeout=timeout, time_precision=time_precision, profile=profile)

if __name__ == '__main__':
    execute()
