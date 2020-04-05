#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Varieties
Litt/Perry Group
Columbia University Mathematics REU 2018
Navtej Singh <singhnav@umich.edu>
'''

# Imports #

import numpy as np

from abc                             import ABCMeta, abstractmethod
from copy                            import deepcopy
from sage.all                        import *
from utils.rational_function         import *
from utils.supersingular             import *
from utils.characters                import *
from utils.bin_lists                 import *
from utils.affine_patch_calculations import *
from utils.gauss_sum_factorization   import *
from utils.reductions                import *
from utils.alpha                     import *

# Classes #

class DefiningEquation():
    def __init__(self, coeffs, powers):
        self.coefficients = coeffs
        self.powers = powers
        self.n_variables = len(coeffs)
        self.degree = lcm(powers)
        if len(coeffs) != len(powers):
            print "ERROR: Ill-defined Equation"

    def to_str(self):
        coeffs = self.coefficients
        powers = self.powers
        if coeffs[0] == 0:
            string = ""
        elif coeffs[0] == 1:
            string = "x_{}^{}".format(0, powers[0])
        else:
            string = "{} * x_{}^{}".format(coeffs[0], 0, powers[0])
        for i in range(1, self.n_variables):
            if coeffs[i] == 0:
                continue
            elif coeffs[i] == 1:
                string += " + x_{}^{}".format(i, powers[i])
            else:
                string += " + {} * x_{}^{}".format(coeffs[i], i, powers[i])
        return string

    def eq_print(self):
        print self.to_str

class AbstractVariety(object):
    '''
    AbstractVariety
    ===============

    Inheritance
    -----------
    ABCMeta, object -> (AbstractVariety) -> ConcreteVariety -> -> Affine Variety, Projective Variety
    '''
    __metaclass__ = ABCMeta

    def __init__(self, p, defining_equation):
        '''
        Inputs
        ------
        <defining_equation> :: list
            Required; the defining equation must be an ordered list of 2-tuples where the
            first element is the coefficient and the second is the exponent.
        '''
        self.characteristic = p
        self.defining_equation = defining_equation
        self.r = self.defining_equation.n_variables - 1

    def print_V(self):
        print "Not Finished"

    def reduce_exponents(self):
        self.defining_equation.powers = exp_reduce(self.powers)

    def resolve_singularities(self):
        self.defining_equation.powers = exp_resolve(self.powers, self.characteristic)

    def reduce_to_base(self):
        self.defining_equation.powers = exp_base(self.powers, self.characteristic)


    @property
    def coeffs(self):
        return self.defining_equation.coefficients

    @property
    def powers(self):
        return self.defining_equation.powers

    @property
    def n_variables(self):
        return self.defining_equation.n_variables

    @property
    def degree(self):
        return self.defining_equation.degree

    def is_supersingular(self):
        return not exist_exceptional_alpha(self.characteristic, self.r, exp_base(self.powers, self.characteristic))

    def exceptional_zeta_factor(self):
        return calc_exceptional_poly(self.characteristic, self.r, self.coeffs, exp_base(self.powers, self.characteristic))

    def newton_polygon(self):
        return compute_newton_polygon(self.r, self.characteristic, exp_base(self.powers, self.characteristic))

    @abstractmethod
    def field(self):
        '''
        N.B. Abstractmethod, must be defined by descendants
        '''
        pass

    @abstractmethod
    def is_smooth(self, k=1):
        '''
        N.B. Abstractmethod, must be defined by descendants
        '''
        pass

class ConcreteVariety(AbstractVariety):
    '''
    ConcreteVariety
    ===============

    Inheritance
    -----------
    AbstractVariety -> (ConcreteVariety) -> -> Affine Variety, Projective Variety
    '''
    def __init__(self, order, defining_equation=None, assigned_field=None):
        '''
        Inputs
        ------
        <order> :: int || Integer
            Required; ConcreteVariety associates a Sage Finite Field object of the given order.
        <defining_equation> :: list
            Optional; the defining equation must be an ordered list of 2-tuples where the
            first element is the coefficient and the second is the exponent.
        <assigned_field> :: Sage Field
            [Future Support] Optional; overrides Finite Field to associate any Sage Field.
        '''
        self.order = order
        factor_list = list(factor(order))
        if len(factor_list) > 1:
            print "ERROR: No field of order {}".format(order)
        else:
            p = factor_list[0][0]
            self.f = factor_list[0][1]
            if defining_equation is not None:
                super(ConcreteVariety, self).__init__(p, defining_equation)
            if assigned_field is not None:
                self.assigned_field = assigned_field

    @property
    def field(self):
        '''
        N.B. Accessible as class attribute 'field'
        '''
        return GF(self.order)

    def set_defining_equation(self, defining_equation):
        super(ConcreteVariety, self).__init__(defining_equation)

    def set_order(self, order):
        self.order = order
        factor_list = list(factor(order))
        if len(factor_list) > 1:
            print "ERROR: No field of order {}".format(order)
        else:
            self.characteristic = factor_list[0][0]
            self.f = factor_list[0][1]

    @abstractmethod
    def num_points(self):
        '''
        N.B. Abstractmethod, must be defined by descendants
        '''
        pass

    def interesting_degree(self):
        if self.degree % self.characteristic == 0:
            return 0
        else:
            return maximal_alpha_num(self.r, self.powers)

    @abstractmethod
    def zeta_degree(self):
        '''
        N.B. Abstractmethod, must be defined by descendants
        '''
        pass

    @abstractmethod
    def picard_number(self):
        '''
        N.B. Abstractmethod, must be defined by descendants
        '''
        pass

    def find_X(self, u):
        r = u.sum()-1
        if sum(u) == 0:
            return 1
        new_coeffs = [self.coeffs[i] for i in range(self.n_variables) if u.is_one(i)]
        new_powers = [self.powers[i] for i in range(self.n_variables) if u.is_one(i)]
        return patch_sols(r, self.characteristic, self.f, new_coeffs, new_powers)

class AffineVariety(ConcreteVariety):
    '''
    AffineVariety
    =============

    Inheritance
    -----------
    ConcreteVariety -> (Affine Variety)
    '''
    def __init__(self, order=None, defining_equation=None):
        '''
        Inputs
        ------
        <order> :: int || Integer
            Optional; if given, ancestor ConcreteVariety associates a Sage Finite Field object of the given order.
        <defining_equation> :: list
            Optional; if given, ancestor AbstractVariety associates the defining equation.
            The defining equation must be an ordered list of 2-tuples where the
            first element is the coefficient and the second is the exponent.
        '''
        if defining_equation is not None and order is not None:
            super(AffineVariety, self).__init__(order, defining_equation)

    def print_V(self):
        print 'Instantiated Affine Variety defined by {} with degree {} and r = {}\nof Order {} (over {}) with Characteristic {} and power {}\n'.format(self.defining_equation.to_str(),
            self.degree, self.r, self.order, self.field, self.characteristic, self.f)

    # Add Affine Variety Methods Here

    def is_smooth(self, k=1):
        for i in self.powers:
            if gcd(self.characteristic, i) > 1:
                return False
        return True

    def num_points(self):
        return patch_sols(self.r, self.characteristic, self.f, self.coeffs, self.powers)

    def zeta_degree(self):
        size = self.interesting_degree()
        num = size
        denom = size + 1
        return num, denom

    def picard_number(self):
        return 2 * count_exceptional_alpha(self.characteristic, self.r, exp_base(self.powers, self.characteristic))[1] + 1

class ProjectiveVariety(ConcreteVariety):
    '''
    ProjectiveVariety
    =================

    Inheritance
    -----------
    ConcreteVariety -> (Projective Variety)
    '''
    def __init__(self, order=None, defining_equation=None):
        '''
        Inputs
        ------
        <order> :: int || Integer
            Optional; if given, ancestor ConcreteVariety associates a Sage Finite Field object of the given order.
        <defining_equation> :: list
            Optional; if given, ancestor AbstractVariety associates the defining equation.
            The defining equation must be an ordered list of 2-tuples where the
            first element is the coefficient and the second is the exponent.
        '''
        if defining_equation is not None and order is not None:
            super(ProjectiveVariety, self).__init__(order, defining_equation)

        self.weights = [self.degree/i for i in self.powers]

    def print_V(self):
        print 'Instantiated Projective Variety defined by {} with degree {} and r = {}\nof Order {} (over {}) with Characteristic {} and power {}\n'.format(self.defining_equation.to_str(),
            self.degree, self.r, self.order, self.field, self.characteristic, self.f)

    # Add Projective Variety Methods Here

    def betti(self, i):
        # the Betti numbers are zero below zero and above the dimension of the variety
        if i < 0 or i > 2*(self.r - 1):
            return 0
        # if i is not the dimension of the variety (where the interesting cohomology is)
        elif i != self.r - 1:
            # return the cohomology of projective space
            if i % 2 == 0:
                return 1
            else:
                return 0
        # i = r-1 the dimension of the variety has interesting cohomology
        else:
            hom = self.interesting_degree()
            if i % 2 == 0:
                return hom + 1
            else:
                return hom

    def zeta_degree(self):
        hom = self.interesting_degree()
        if self.r % 2 == 0:
            num = hom
            denom = self.r
        else:
            num = 0
            denom = self.r + hom
        return num, denom

    def picard_number(self):
        return count_exceptional_alpha(self.characteristic, self.r, exp_base(self.powers, self.characteristic))[1] + self.r

    def euler_char(self):
        num, denom = self.zeta_degree()
        return denom - num

    def is_smooth(self, k=1):
        p = self.characteristic
        q = self.order

        for i in self.powers:
            if gcd(i, q) > 1:
                return False
        L = [l for (l, e) in list(factor(lcm(self.weights)))]
        #print w, L
        for l in L:
            if gcd(l, q - 1) == 1:
                continue
            else:
                which_vars = [0]*(self.r + 1)
                for i in range(len(which_vars)):
                    if gcd(l, self.weights[i]) > 1:
                        which_vars[i] = 1
                u = bin_list(which_vars)
                if gcd(l, q - 1) > 1 and self.find_X(u) > 1:
                    return False
        return True

    # check if a weighted projective variety is actually defined in projective space
    # meaning the weights are all 1 or equivalently the defining equation is homogeneous
    def is_projective(self):
        for w in self.weights:
            if w != 1:
                return False
        return True

    def num_points(self):
        if self.is_projective():
            return (patch_sols(self.r, self.characteristic, self.f, self.coeffs, self.powers) - 1)/ (self.order - 1)
        else:
            u = all_one(self.r)
            X = dict()

            for t in u.all_less_than():
                X[t.to_int()] = self.find_X(t)

            C = dict()

            for t in u.all_less_than():
                C[t.to_int()] = alt_sum_over_bins(X, t) + (-1)**(t.sum())
            N = 0

            for t in u.all_less_than():
                N += C[t.to_int()] * Integer(gcd(find_d_bin(t, self.powers), self.order - 1)) / Integer(self.order - 1)

            return N

# Functions #

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

def get_solutions(V, m):
    q = V.order
    V_prime = deepcopy(V)
    sols = []
    for k in range(1, m):
        V_prime.set_order(q**k)
        sols.append(V_prime.num_points())
    return sols

def zeta_function(V):
    num, denom = V.zeta_degree()
    char = denom - num
    sols = get_solutions(V, 2*max(num, denom) + 4)
    return find_rational(sols, 1, max(-char + 1, 0), denom + 1)

def format_zeta_function(V):
    ratl_fn = zeta_function(V)
    numerator = print_rational(ratl_fn[0], V.order)
    denominator = print_rational(ratl_fn[1], V.order)
    return "\n{}\n-------------------\n{}".format(numerator, denominator)

def affine_zeta_prime_subfield(V):
    return zeta_function(AffineVariety(V))

def projective_zeta_prime_subfield(V):
    return zeta_function(ProjectiveVariety(V))

def is_variety_supersingular(V):
    ratl_fn = zeta_function(V)
    return is_supersingular(ratl_fn[0], ratl_fn[1], V.order)

    poly = zeta_function(V)[1]
    if poly == 0:
        print "ERROR: zeta function not found"
        return None
    # in the case of a projective surface we ignore the top and bottom homology
    if isinstance(V, ProjectiveVariety) and V.r == 3:
        R = QQ['T']
        t = R.gen()
        poly = R(poly/(1 - t)/(1 - V.order**2 * t))
    return poly.change_ring(Qp(V.characteristic)).newton_polygon().vertices()

def analyze_variety(V):
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

# Tests

def test_affine(order, defining_equation):
    av = AffineVariety(order, defining_equation)
    av.print_V()
    print "Is this variety smooth? {}\n".format(av.is_smooth())
    print "Number of Solutions: {}\n".format(get_solutions(av, 10))
    print "Zeta Function: {}\n".format(format_zeta_function(av))
    print "Newton Polygon: {}\n".format(newton_polygon(av))

def test_projective(order, defining_equation):
    pv = ProjectiveVariety(order, defining_equation)
    pv.print_V()
    print "Is this variety smooth? {}\n".format(pv.is_smooth())
    print "Number of Solutions: {}\n".format(get_solutions(pv, 10))
    print "Zeta Function: {}\n".format(format_zeta_function(pv))
    print "Newton Polygon: {}\n".format(newton_polygon(pv))

# Main

def main(affine_order=3, affine_defining_equation=DefiningEquation([1,1,1,1], [5,5,5,5]),
    projective_order=3, projective_defining_equation=DefiningEquation([1,1,1,1], [5,5,5,5])):
    test_affine(affine_order, affine_defining_equation)
    test_projective(projective_order, projective_defining_equation)

def tests():
    for p in [2,3,5,7,11,13,17,19,23]:
        for n in range(4, 40):
            V = AffineVariety(p, DefiningEquation([1,1,1,1], [n,n,n,n]))
            if n % p != 0:
                f = Mod(p, n).multiplicative_order()
                if f % 2 == 0:
                    print "p = {}, n = {}, f = {}, p**(f/2) ~ {} mod {} and SS? {}".format(p, n, f, Mod(p**(f/2), n), n, V.is_supersingular())
                else:
                    print "p = {}, n = {}, f = {} and SS? {}".format(p, n, f, V.is_supersingular())

# Execution #

if __name__ == '__main__':
    pass
