#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Rational Function
Litt/Perry Group
Columbia University Mathematics REU 2018
'''

# Imports #

from sage.all import *
from recurrence import *

# Functions #

def zeta_fn(N_sols, m):
    '''
    list_to_zeta
    Input: N_sols, m = number of terms in N_sols
    Output: Zeta polynomial ''
    '''
    R = QQ['T']
    T = R.gen(0)
    poly = Integer(0) * T
    for i in range(1, m + 1):
        poly += Integer(N_sols[i - 1]) / Integer(i) * (T ** i)
    zeta = Integer(1) + Integer(0) * T
    poly_power = poly
    for j in range(1, m + 1):
        zeta += 1 / (Integer(factorial(j))) * poly_power
        zeta = zeta % (T ** (m + 1))
        poly_power *= poly
        poly_power = poly_power % (T**(m+1))
    return zeta

def find_rational(data, check_limit=4, start=0, dim=2):
    '''
    find_rational(recur, zeta_list)
    Input: zeta - as above, m - as above
    Output: rational function for the zeta polynomial if found, (0,0) else
    '''
    m = len(data)
    zeta = zeta_fn(data, m)
    #print(zeta)
    S = QQ['T']
    T = S.gen(0)
    zeta_list =  zeta.list()[:m+1] #First m+1 terms of zeta polynomial
    recur = find_recurrence(zeta_list, check_limit, start, dim)
    if (recur == 0): #If no recurrence was found
        return (0,0)
    recur = recur[::-1] #Flip the recurrence list
    exp_q = len(recur)
    q = -1

    #Loop that finds the denominator (q)
    for i in range(exp_q):
        q += Rational(recur[i]) * (T**(i+1))
    coeffs = zeta.list()[:len(recur) + start + 1]

    #Computing the numerator (p) by multiplying first len(recur) + 1 terms of zeta by q
    zeta_copy = Integer(0) * T
    for index, coeff in enumerate(coeffs):
        zeta_copy += coeff*(T**index)
    p = zeta_copy * q
    p_coeffs = p.list()
    p = Integer(0) * T
    for i in range(len(recur) + start + 1):
        p += p_coeffs[i] * T**i
    return (p, q)

def print_rational(func, prim_power):
    S = QQ['T']
    T = S.gen(0)
    if (func == 0):
        # print("0")
        return 0
    factored = func.factor()
    constant = func.list()[-1]
    print_string = ''
    for exprs in factored:
        poly = exprs[0]
        exp = exprs[1]
        shift = log(abs(poly(0)), prim_power)/poly.degree()
        shift_power = poly(0)**(1/poly.degree())
        prim_exp = floor(shift)
        shift_power = Integer(prim_power)**Integer(prim_exp)
        prim_poly = poly(shift_power * T)/(shift_power**(poly.degree()))
        constant *= shift_power**(exp * poly.degree())
        print_string += '['
        for index, coeff in enumerate(prim_poly.list()):
            if (index == 0):
                print_string += '%s'%(coeff)
            if (coeff < 0 and index != 0):
                print_string += '%s*'%(coeff)
            if (coeff > 0 and index != 0 and coeff != 1):
                print_string += '+%s*'%(coeff)
            if (coeff > 0 and index != 0 and coeff == 1):
                print_string += '+'
            if (coeff != 0 and index == 1 and prim_exp != 0):
                print_string += '%s^(%s)*T'%(prim_power, -prim_exp) if -prim_exp != 1 else '%s*T' % prim_power
            if (coeff != 0 and index == 1 and prim_exp == 0):
                print_string += 'T'
            if (coeff != 0 and index > 1 and prim_exp != 0):
                if -prim_exp != 1 and index != 1:
                    print_string += '(%s^(%s)*T)^%s' % (prim_power, -prim_exp, index)
                elif -prim_exp != 1 and index == 1:
                    print_string += '(%s^(%s)*T)' % (prim_power, -prim_exp)
                elif -prim_exp == 1 and index != 1:
                    print_string += '(%s*T)^%s' % (prim_power, index)
                elif -prim_exp == 1 and index == 1:
                    print_string += '(%s*T)' % prim_power
            if (coeff != 0 and index > 1 and prim_exp == 0):
                print_string += 'T^%s'%(index)
        print_string += ']^%s'%(exp) if exp > 1 else ']'
    if len(print_string) > 0:
        if constant > 1:
            return '%s*%s' % (constant , print_string)
        else:
            return '%s' % print_string
    else:
        return '%s' % constant

def main():
    data = [14641, 214358881, 3138428376721, 45949729863572161, 672749994932560009201, 9849732675807611094711841,
    144209936106499234037676064081, 2111377674535255285545615254209921, 30912680532870672635673352936887453361,
    452592555681759518058893560348969204658401, 6626407607736641103900260617069258125403649041]
    result = find_rational(data)
    print '{}, {}\n'.format(result[0], result[1])

# Execution #

if __name__ == '__main__':
    main()
