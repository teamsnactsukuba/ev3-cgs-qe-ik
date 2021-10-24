from math import sqrt
from ctypes import * 
from sympy import *
from sympy import im
from sympy import polys
from scipy.optimize import fsolve
from fractions import Fraction
from decimal import Decimal
import random
import time
import numpy as np
import sympy as sp

SamplePoint = []

# Although (X,Y,Z) should satisfy 
# X:(-310,310), Y:(-310,310), Z:(0,414),
# in the actural computation, many points are detected as infeasible;
# thus, take (X,Y,Z) from smaller region
# The condition of if sentence is brought from the distance from Joint 4
for J in range(3):
    X = Fraction(random.uniform(-200,200))
    Y = Fraction(random.uniform(-200,200))
    Z = Fraction(random.uniform(0,300))
    X = X.limit_denominator(100)
    Y = Y.limit_denominator(100)
    Z = Z.limit_denominator(100)
    if (X-(44*sqrt(2))/(sqrt(X**2+Y**2)))**2+(Y-(44*sqrt(2))/(sqrt(X**2+Y**2)))**2+(Z-104-44*sqrt(2))**2 <= (248)**2 and X**2+Y**2 != 0 and not [X,Y,Z] in SamplePoint :
        SamplePoint.append([X,Y,Z])


print(SamplePoint)
print(len(SamplePoint))