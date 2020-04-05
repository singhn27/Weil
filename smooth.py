# Imports #

from sage.all          import *
from utils.bin_lists   import *
from utils.characters        import *
from itertools         import product
from point_counting    import find_X

def is_diagonal(monomial):
    found = False
    for i in monomial:
        if found and i != 0:
            return False
        elif not found and i != 0:
            found = True
    return True

def all_monomials(powers):
    d = lcm(powers)
    w = [d/i for i in powers]
    r = len(powers)-1
    lists = [range(n+1) for n in powers]
    possible = product(*lists)
    monomials = []
    for a in possible:
        if sum([x*y for x,y in zip(a, w)]) == d and not is_diagonal(a):
            monomials.append(a)
    return monomials

def find_constraint(powers, monomial):
    b = lcm([n/gcd(n,a) for (n,a) in zip(powers, monomial)])
    c = [a*b/n for (n,a) in zip(powers, monomial)]
    return (b,c)


def WP_smooth(powers, q):
    r = len(powers)-1
    monomials = all_monomials(powers)
    N = len(monomials) + r
    string = "x0"
    for j in range(1, r+1):
        string = string + " + x{}".format(j)

    constrain_eqs = [string]

    for i, m in enumerate(monomials):
        b,c = find_constraint(powers, m)
        string = "x0**{}".format(c[0])
        for j in range(1, r+1):
            string = string + " * x{}**{}".format(j, c[j])
        constrain_eqs.append(string + " - x{}^{}".format(i+r+1, b))
        print m, string + " - x{}^{}".format(i+r+1, b)
    Proj = ProjectiveSpace(N, GF(q))
    print Proj
    print constrain_eqs
    V = Proj.subscheme(constrain_eqs)
    print V
    return V.is_smooth()

def my_WP_smooth(coeffs, powers, p, q, k=1):
    r = len(powers)-1
    d = lcm(powers)
    w = [d/i for i in powers]
    f = log(q, p)

    for i in powers:
        if gcd(i,q) > 1:
            return False
    L = [l for (l, e) in list(factor(lcm(w)))]
    #print w, L
    for l in L:
        if gcd(l, q-1) == 1:
            continue
        else:
            which_vars = [0]*(r+1)
            for i in range(len(which_vars)):
                if gcd(l, w[i]) > 1:
                    which_vars[i] = 1
            u = bin_list(which_vars)
            #u.bin_print()
            gauss_sums = get_sums(p, f, k, powers, r) # Compute 2-dim matrix of gauss sums
            chars_of_coeffs = get_chars(p, f, k, powers, coeffs, r) # Compute 2-dim matrix of chi_{alpha_i}(a_i^{-1})
            if gcd(l, q-1) > 1 and find_X(u, p, q, 1, coeffs, powers, gauss_sums, chars_of_coeffs, k)[0] > 1:
                return False
    return True

if __name__ == '__main__':
    coeffs = [1,1,1,1]
    powers = [4,4,8,8]
    p = 13
    q = 13

    def in_order(ps):
        for i in range(len(ps)-1):
            if ps[i] > ps[i+1]:
                return False
        return True

    lists = [range(1, 11) for n in range(0,4)]
    possible = product(*lists)
    for v in possible:
        if in_order(v):
            for p in [2,3,5,7,11,13,17,19,23]:
                print "{}  over p = {} is smooth? {}".format(v, p, my_WP_smooth(coeffs, v, p, p, 1))
