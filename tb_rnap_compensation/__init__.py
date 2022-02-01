#! /usr/bin/env python3

import numpy

def calculate_fisher_pvalue(array):

    assert isinstance(array, numpy.array), "have to pass a numpy array!"

    assert array.shape = (2,2), "has to be a 2x2 array!"

    # and use fisher.pvalue
    p = pvalue(array[0,0],array[0,1],array[1,0],array[1,1])

    return(p)
