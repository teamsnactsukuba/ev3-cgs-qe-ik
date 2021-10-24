#!/usr/bin/env python
# coding: UTF-8
from sympy import *

t1 = Symbol('t1')
t4 = Symbol('t4')
t7 = Symbol('t7')
x = Symbol('x')
s1 = Symbol('s1')
s2 = Symbol('s2')
s3 = Symbol('s3')
c1 = Symbol('c1')
c2 = Symbol('c2')
c3 = Symbol('c3')

def DHMatrix(a,A,d,T):
        return Matrix([
                [cos(T),-sin(T),0,a],
                [cos(A)*sin(T),cos(A)*cos(T),-sin(A),-d*sin(A)],
                [sin(A)*sin(T),sin(A)*cos(T),cos(A),d*cos(A)],
                [0,0,0,1]
        ])
        
M01 = DHMatrix(0,0,80,t1)
M12 = DHMatrix(0,0.5*pi,0,0.25*pi)
M23 = DHMatrix(88,0,0,pi/4.0)
M34 = DHMatrix(24,0,0,t4)
M45 = DHMatrix(96,0,0,-pi/2.0)
M56 = DHMatrix(16,0,0,pi/2.0)
M67 = DHMatrix(40,0,0,t7)
M78 = DHMatrix(112,0,0,0)
    
M = M01*M12*M23*M34*M45*M56*M67*M78
f1 = M[0,3]
f2 = M[1,3]
f3 = M[2,3]

