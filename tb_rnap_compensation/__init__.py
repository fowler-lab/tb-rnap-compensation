#! /usr/bin/env python3

import numpy

from fisher import pvalue
from scipy.stats import chi2_contingency

import decimal
def calculate_fisher_pvalue(array):

    assert isinstance(array, numpy.ndarray), "have to pass a numpy array!"

    assert array.shape == (2,2), "has to be a 2x2 array!"

    # and use fisher.pvalue
    p = pvalue(array[0,0],array[0,1],array[1,0],array[1,1])

    return(p)

def calculate_chi_square_pvalue(array):

    assert isinstance(array, numpy.ndarray), "have to pass a numpy array!"

    assert array.shape == (2,2), "has to be a 2x2 array!"

    stat, p, dof, expected = chi2_contingency(array)

    return(stat, p, dof, expected)

def numerical_test(n = 10000, N = 10000, res_obs = 14000, other_obs = 40000, both_obs = 8800):
    
    p_res = res_obs/ 64722
    p_other = other_obs/ 64722
    n_both_obs = (both_obs/ 64722) * N
    
    n_both = []
    
    for i in range(n):

        N = N

        resistant = numpy.random.choice(a = [True, False], size = N, p = [p_res, 1-p_res])

        other = numpy.random.choice(a = [True, False], size = N, p = [p_other, 1-p_other])

        both = other[resistant]

        n_both.append(both.sum())

    p_value = decimal.Decimal((numpy.array([n_both]) >= n_both_obs).sum().item())/ decimal.Decimal(n)
    
    return p_value