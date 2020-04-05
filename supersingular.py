min_phis = [0, 2, 6, 6, 12, 12, 18, 18, 30, 30, 30, 30, 42, 42, 42, 42, 60, 60, 60, 60, 66, 66, 66, 66, 90, 90, 90, 90, 90, 90, 90, 90, 120, 120, 120, 120, 126, 126, 126, 126, 150, 150, 150, 150, 150, 150, 150, 150, 210, 210, 210, 210, 210, 210, 210, 210, 210, 210, 210, 210, 210, 210, 210, 210, 240, 240, 240, 240, 240, 240, 240, 240, 270, 270, 270, 270, 270, 270, 270, 270, 330, 330, 330, 330, 330, 330, 330, 330, 330, 330, 330, 330, 330, 330, 330, 330, 420, 420, 420, 420] #min_phis[i] = minimal n such that phi(N) > i for all N > n 

from sage.all import *

def cyclo_ub(i):
	'''
	Input: The degree of a polynomial over QQ
	Output: Finds an upper bound for the degree of the largest cyclotomic poly dividing 
	'''
	if i < 100:
		return min_phis[i]
	return floor(i * log(i))

def is_factor_super(fact, q):
	'''
	Input: An irreducible factor (fact) over QQ, the order of the finite field q
	Output: 1 if all the roots are of the form q^(i/2) * root of unity, 0 else 
	'''
	R = QQ['T']
	T = R.gen(0)
	root_mag = fact(0)**(1/fact.degree())
	scaled_fact =1/fact(0) * fact(root_mag * T)
	for i in range(1, cyclo_ub(fact.degree()) + 1):
		if gcd(T**i - 1, scaled_fact) != 1:
			return 1
	return 0
	

def is_supersingular(num, den, q):
	'''
	Input: The numerator (num) and denominator (den) of a rational function, finite field order q
	Output: True if num/den is supersingular over F_q, False else
	'''
	if num == 0 or den == 0:
		return None
	num_factor = num.factor()
	den_factor = den.factor()
	for fact_pair in num_factor:
		if not is_factor_super(fact_pair[0],q):
			return False
	for fact_pair in den_factor:
		if not is_factor_super(fact_pair[0],q):
			return False
	return True

