# -*- coding: utf-8 -*-

import sympy as sp

# Used for real root counting with the Descartes' rule of signs

def descaltes4(H):

    # Count the number of sign changes of the nonzero elements in sequence
    # Input: H = [h1, ... , hn]

    F = [] # excluding zeros in H
    for P in H:
        a = sp.sign(P)
        if a != 0:
            F.append(a)
    I = 1
    C = 0
    M = len(F)
    N = M
    while N > 1:
        if F[I-1] == 1 and F[I] == -1:
            C += 1
        if F[I-1] == -1 and F[I] == 1:
            C += 1
        I += 1
        N -= 1
    N = M

    F2 = [] 
    # for F = [f1, ..., fm], calculate F2 as
    # F2 = [f1, -f2, f3, -f4, ... , (-1)^m fm]
    for J in range(N):
        if J % 2 ==0:
            F2.append(F[J])
        else :
            F2.append(F[J]*(-1))
    I = 1
    D = 0
    while N > 1:
        if F2[I-1] == 1 and F2[I] == -1:
            D += 1
        if F2[I-1] == -1 and F2[I] == 1:
            D += 1
        I += 1
        N -= 1
    return C - D
