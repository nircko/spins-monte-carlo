import numpy as np
from math import pi, cos, sqrt

def coth(x):
    y=(np.cosh(x) / np.sinh(x))
    return y


def brillouin(j,x):
    return (1 / j) * ((j+1/2) * (coth((j+1/2) * x) - (1/2) * coth(x/2)))

def chi_anal(x):
    return  0.25*(1/np.sinh(x/2))**2-(1/np.sinh(x))**2



#import sys,os, inspect

#sys.path.append(r'C:\Users\Nir Goldfriend\PycharmProjects\statistical_numeric_ex2017')
#import moduleExNum
