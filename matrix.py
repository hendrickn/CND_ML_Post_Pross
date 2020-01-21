from __future__ import absolute_import, division, print_function, unicode_literals
from numba import vectorize, cuda
import math


#import tensorflow as tf

@vectorize(['float64(float64, float64)'],nopython=True, target='cpu')
def dot_mul(a,b) :
    return a*b


@vectorize(['float64(float64)'],nopython=True, target='cpu')
def dot_exp(a) :
    return math.exp(a)



