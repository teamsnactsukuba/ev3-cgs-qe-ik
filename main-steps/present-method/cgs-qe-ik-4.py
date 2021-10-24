#! /usr/bin/env python3

# cgs-qe-ik-0.py
# Inverse kinematic computation 
# with Quantifier Elimination based on Comprehensive Groebner Systems (CGS-QE)
# for solving inverse kinematic problems
# base varsion

# Copyright (C) Team SNAC Tsukuba, 2021

from math import sqrt, sin, cos
from ctypes import * 
from fractions import Fraction
import time
import numpy as np
import sympy as sp
import sys

import platform
import os
import datetime

import cgs
import partition
import descartes
import SamplePoint

print (datetime.datetime.now())

# SymPy initializations
# (x, y, z): coordinates of the end effector
# s_i = sin(theta_i)
# c_i = cos(theta_i)
# t: the variable in the characteristic polynomial of the Hermite quadratic form

x = sp.Symbol('x')
y = sp.Symbol('y')
z = sp.Symbol('z')
s_1 = sp.Symbol('s_1')
c_1 = sp.Symbol('c_1')
s_4 = sp.Symbol('s_4')
s_7 = sp.Symbol('s_7')
c_4 = sp.Symbol('c_4')
c_7 = sp.Symbol('c_7')
t = sp.Symbol('t')

# A: a list of partitions in the CGS
# B: a list of Groebner bases of the CGS
# E: a list of degree of each polyonials in the CGS
# Points: a list of sample points
A = partition.A
B = cgs.B
E = cgs.E
Points = SamplePoint.Points

# list_CellNum: a list of index of partition in A which contains the sample point
# list_CellNumTime: a list of computing time for calculating the index of partition in A which contains the sample point
# list_CountRealRoots: a list of the number of real roots
# list_CountRealRootsTime: a list of computing time for counting the number of real roots
# list_SubsTime: a list of computing time for substituting variables with the sample point in counting the number of real roots
# list_DescarteTime: a list of computing time for the Desrates' rule of signs 
# list_SolveTime: a list of computing time for solving a system of algebraic equations
# list_IKSolvingTime: a list of computing time for solving inverse kinematic problem
# list_Error: a list of error of the solution of the inverse kinematic computation from the given sample point
# list_ErrorTime: a list of computing time for estimating errors

list_CellNum = []
list_CellNumTime = []
list_CountRealRoots = []
list_CountRealRootsTime = []
list_SubsTime = []
list_DescarteTime = []
list_ValidSamplePoint = []
list_SolveTime = []
list_IKSolvingTime = []
list_Error = []
list_ErrorTime = []

# OpenXM initializations

# An environment variable OpenXM_HOME must be defined

OpenXM_HOME = os.environ['OpenXM_HOME']    
print("OpenXM_Home =", OpenXM_HOME) 
pf = platform.system()
print("Platform =", pf)
if pf == 'Linux':
    LIBGMP = OpenXM_HOME + "lib/libgmp.so"   
    LIBGC = OpenXM_HOME + "lib/libgc.so"     
    LIBMPFR = OpenXM_HOME + "lib/libmpfr.so" 
    LIBOX = OpenXM_HOME + "lib/libox.so"     
elif pf == 'Darwin':
    LIBGMP = OpenXM_HOME + "lib/libgmp.dylib"
    LIBGC = OpenXM_HOME + "lib/libgc.dylib"
    LIBMPFR = OpenXM_HOME + "lib/libmpfr.dylib"
    LIBOX = OpenXM_HOME + "lib/libox.dylib"
    
OX_SERVER_HOST = create_string_buffer("127.0.0.1".encode('utf-8'))
OX_PROG1 = create_string_buffer("ox".encode('utf-8'))
OX_PROG2 = create_string_buffer("ox_asir".encode('utf-8'))    

libox=CDLL("libgmp.so", mode=RTLD_GLOBAL)
libox=CDLL("libgc.so", mode=RTLD_GLOBAL)
libox=CDLL("libmpfr.so", mode=RTLD_GLOBAL)
libox=CDLL("libox.so", mode=RTLD_GLOBAL)

# from OpenXM_HOME/src/ox_toolkit/ox_toolkit.h
# typedef struct {
#     int tag;
#     int flag;
# } table;
class TABLE (Structure):
    _fields_ = [
        ("tag", c_int),
        ("flag", c_int)
        ]

# from OpenXM_HOME/src/ox_toolkit/ox_toolkit.h
# typedef struct mathcap {
#     table *cmotbl;
#     table *smtbl;
#     char  **opts;
# } mathcap;
class MATHCAP (Structure):
    _fields_ = [
        ("cmotbl", POINTER(TABLE)),
        ("smtbl", POINTER(TABLE)),
        ("opts", POINTER(c_char_p))
        ]
        
# from OpenXM_HOME/src/ox_toolkit/ox_toolkit.h
# typedef struct OXFILE{
#     int fd;
#     int (*send_int32)(struct OXFILE *oxfp, int int32);
#     int (*receive_int32)(struct OXFILE *oxfp);
#     int serial_number;
#     int received_serial_number;
#     struct OXFILE *control;  /* pointer to his control server. */
#     struct mathcap *mathcap;
#     int error;
#     char *wbuf;
#     int wbuf_size;
#     int wbuf_count;
#     int (*send_double)(struct OXFILE *oxfp, double int64);
#     double (*receive_double)(struct OXFILE *oxfp);
# } OXFILE;
# Since OXFILE has self referential structure, define the class with 'incomplete types'
# https://docs.python.org/ja/3/library/ctypes.html#incomplete-types
class OXFILE (Structure):
    pass

libox.send_int32.restype = POINTER(c_int)
libox.send_int32.argtypes = (POINTER(OXFILE), c_int)
libox.receive_int32.restype = c_int
libox.receive_int32.argtype = POINTER(OXFILE)
libox.send_double.restype = c_int
libox.send_double.argtypes = (POINTER(OXFILE), c_double)
libox.receive_double.restype = c_double
libox.receive_double.argtype = POINTER(OXFILE)

OXFILE._fields_ = [
    ("fd", c_int),
    ("send_int32", POINTER(c_int)),
    ("receive_int32", POINTER(c_int)),
    ("serial_number", c_int),
    ("received_serial_number", c_int),
    ("control", POINTER(OXFILE)),
    ("mathcap", POINTER(MATHCAP)),
    ("error", c_int),
    ("wbuf", POINTER(c_char)),
    ("wbuf_size", c_int),
    ("wbuf_count", c_int),
    ("send_double", POINTER(c_int)),
    ("receive_double", POINTER(c_double))
    ]
    
libox.ox_start.restype = POINTER(OXFILE)
libox.ox_start.argtypes = (POINTER(c_char), POINTER(c_char), POINTER(c_char))
libox.ox_reset.restype = None
libox.ox_reset.argtype = POINTER(OXFILE)
libox.ox_shutdown.restype = None
libox.ox_shutdown.argtype = POINTER(OXFILE)
libox.ox_execute_string.restype = None
libox.ox_execute_string.argtypes = (POINTER(OXFILE), c_char_p)
libox.ox_popString.restype = c_char_p
libox.ox_popString.argtype = POINTER(OXFILE)

# Obtain current director for OpenXM

pwd = os.getcwd()

# Start the ox server

sv5 = libox.ox_start(OX_SERVER_HOST, OX_PROG1, OX_PROG2)

# Asir: load the program and the list of characteristic polynomials

C_DAT_FILE = "\"" + pwd + "/C.dat" + "\""
ox_command_string = "C=bload(" + C_DAT_FILE + ")$"
print("ox_command_string:", ox_command_string)
OX_COMMAND = create_string_buffer(ox_command_string.encode('utf-8'))
libox.ox_execute_string(sv5, OX_COMMAND)

SUBST_CHAR_POLY_FILE = "\"" + pwd + "/subst_char_poly.rr" + "\""
ox_command_string = "load(" + SUBST_CHAR_POLY_FILE + ")$"
print("ox_command_string:", ox_command_string)
OX_COMMAND = create_string_buffer(ox_command_string.encode('utf-8'))
libox.ox_execute_string(sv5, OX_COMMAND)

r = range(300, 400) # range of sample points

for k in r:

    print("-----")
    print("Case:", k)

    # Take a point from SamplePoint.py
    # Let (x, y, z) = (e1, e2, e3) be each coordinates

    Point = SamplePoint.Points[k]
    e1 = Point[0]
    e2 = Point[1]
    e3 = Point[2]

    print("X =", e1)
    print("Y =", e2)
    print("Z =", e3)

    TimeStart = time.time()
    
    # Choose a partition in A that containts the sample point (e1, e2, e3)
    # by testing the following conditions:
    # for a partition expressed as X = [W, V]
    # where W = [w1, ... , wk] and V = [v1, ... , vl] with w1, ..., wk, v1, ... , vl are polynomials in x, y, z,
    # w1 = ... wk = 0 and v1 != 0, ..., vl != 0
    # by substituting x=e1, y=e2, z=e3
    
    partition_index = -1
    num = 0
    bre = 0
    for X in A:
        zero_list0 =[0] * len(X[0])
        cell_list0 =[]
        for W in X[0]:
            if W != 0:
                gen_0 = W.subs([(x,e1),(y,e2),(z,e3)])
                cell_list0.append(gen_0)
            else :
                cell_list0.append(W)
        if cell_list0 == zero_list0:
            zero_list1 = [0] * len(X[1])
            cell_list1 = []
            for V in X[1]:
                if V != 1:
                    gen_1 = V.subs([(x,e1),(y,e2),(z,e3)])
                    cell_list1.append(gen_1)
                else :
                    cell_list1.append(V)
            if cell_list1 != zero_list1:
                partition_index = num
                bre = 1
                break
        if bre == 1:
            break
        num += 1

    # Put the index of appropriate partition into list_CellNum
    print("The index of partition:", partition_index)
    TimeSelectPartition = time.time()
    list_CellNum.append(partition_index)
    list_CellNumTime.append(TimeSelectPartition - TimeStart)

    # Choose the corresponding characteristic polynomial of the Hermite quadratic form,
    # for counting the number of real roots

    ox_command_string = "V = subst_char_poly_vect (C, " + str(partition_index) + ", " + str(e1) + ", " + str(e2) + ", " + str(e3) + ")$"
    print("ox_command_string:", ox_command_string)
    OX_COMMAND = create_string_buffer(ox_command_string.encode('utf-8'))
    libox.ox_execute_string(sv5, OX_COMMAND)

    # Extract the size of the vector of the coefficients in the characteristic polynomial
    # for V = [c_m, ... , c_0] , let chi2_length := m

    ox_command_string = "length(V);"
    OX_COMMAND = create_string_buffer(ox_command_string.encode('utf-8'))
    libox.ox_execute_string(sv5, OX_COMMAND)
    OX_OUTPUT = libox.ox_popString(sv5)
    chi2_length = int(OX_OUTPUT)

    chi2 = []

    # Substitute x=e1, y=e2, z=e3
    
    # for V in chi:
    #     if V !=1 :
    #         chi2.append(V.subs([(x,e1),(y,e2),(z,e3)]))
    #     else :
    #         chi2.append(1)

    for i in range(chi2_length):
        ox_command_string = "V[" + str(i) + "];"
        OX_COMMAND = create_string_buffer(ox_command_string.encode('utf-8'))
        libox.ox_execute_string(sv5, OX_COMMAND)
        OX_OUTPUT = libox.ox_popString(sv5)
        f = OX_OUTPUT.decode("utf-8")
        f = f.replace("^","**")
        f = f.replace("(2)^(1/2)", "sqrt(2)")
        # chi2.append(sympify(f))
        chi2.append(f)

    TimeSubstitution = time.time()
    
    # Count the number of real roots with the Descarts' rule of signs
    
    boo = descartes.descaltes4(chi2)
    list_CountRealRoots.append(boo)
    TimeDescartes = time.time()
    
    print ("Time for real roots counting:", TimeDescartes - TimeSelectPartition)

    list_SubsTime.append(TimeSubstitution - TimeSelectPartition)
    list_DescarteTime.append(TimeDescartes - TimeSubstitution)
    list_CountRealRootsTime.append(TimeDescartes - TimeSelectPartition)

    # If there are no real roots, stop the computation for this sample points here

    if boo == 0 :
        print("no real root!")
        break

    # If there exist real roots, then choose the Groebner basis corresponding the partition

    list_ValidSamplePoint.append(k)
    gb = B[partition_index]
    deg_g = E[partition_index]

    # The number of variables change for the values of (e1, e2) = (x,y):
    # e1 != 0 or e2 !=0: (c_1,s_1,c_4,s_4,c_7,s_7)
    # e1 == 0 and e2 == 0: (c_4,s_4,c_7,s_7)
    list_Solves = []


    if e1 == 0 and e2 == 0:

        # In the case e1 == e2 == 0 (corresponding to s_1 == 0 and c_1 == 1),
        # use g0, g1, g2, g3 for solving the inverse kinematic problem
        # for deriving s_7, c_7, s_4, c_4 

        g0 = gb[0].subs(z,e3)
        g1 = gb[1].subs(z,e3)
        g2 = gb[2].subs(z,e3)
        g3 = gb[3].subs(z,e3)
        lc_g1 = g1.coeff(c_7,deg_g[1])
        lc_g2 = g2.coeff(s_4,deg_g[2])
        lc_g3 = g3.coeff(c_4,deg_g[3])

        deg_g0 = deg_g[0]
        coef_0 = []

        # Solve g0 for s_7

        for W in range(deg_g0+1):
            n = g0.coeff(s_7,W)
            n = n.evalf(20)
            coef_0.insert(0,n)
    
        list_s7 = np.roots(coef_0)
        num_s7 = len(list_s7)
        # Eliminate imaginary roots
        fake_array = []
        for P in range(num_s7):
            if sp.im(list_s7[P]) ==0:
                fake_array.append(list_s7[P])
        list_s7 = fake_array
        num_s7 = len(list_s7)

        # Eliminate multiplicity from multiple roots
        # list -> dict -> list conversion 
        # e.g. [1,1,2]--->[1,2]
        fake_array2 = list(dict.fromkeys(list_s7))
        list_s7 = fake_array2
        num_s7 = len(list_s7)

        # Verify to apply the extension theorem
        fake_array3 = []
        for P in range(num_s7):
            if  lc_g1.subs(s_7,list_s7[P]) == 0:
                print("lc(g1)=0!")
            else:
                fake_array3.append(list_s7[P])
        list_s7 = fake_array3
        num_s7 = len(list_s7)


        for Q in range(num_s7):
            # Put the root s_7 into g1; then, g1 becomes a univariate polynomial in c_7
            # Solve g1 for c_7
            g1_subs = g1.subs(s_7, list_s7[Q])
            deg_g1 = deg_g[1]
            coef_1 = []
            for P in range(deg_g1 + 1):
                n = g1_subs.coeff(c_7, P)
                n = n.evalf(20)
                coef_1.insert(0,n)
            list_c7 = np.roots(coef_1)
            num_c7 = len(list_c7)

            fake_array = []
            for P in range(num_c7):
                if sp.im(list_c7[P]) ==0:
                    fake_array.append(list_c7[P])
            list_c7 = fake_array
            num_c7 = len(list_c7)
        
            fake_array2 = list(dict.fromkeys(list_c7))
            list_c7 = fake_array2
            num_c7 = len(list_c7)
            
            fake_array3 = []
        
            for P in range(num_c7):
                if  lc_g2.subs(c_7,list_c7[P]) == 0:
                    print("lc(g2)=0!")
                else:
                    fake_array3.append(list_c7[P])
        
            list_c7 = fake_array3
            num_c7 = len(list_c7)

            for R in range(num_c7):
                # Put s_7 and c_7 into g2; then, g2 becomes a univariate polynomial in s_4
                # Solve g2 for s_4
                g2_subs = g2.subs([(s_7,list_s7[Q]),(c_7,list_c7[R])])
                deg_g2 = deg_g[2]
                coef_2 = []
                for P in range(deg_g2 + 1):
                    n = g2_subs.coeff(s_4, P)
                    n = n.evalf(20)
                    coef_2.insert(0,n)
                list_s4 = np.roots(coef_2)
                num_s4 = len(list_s4)

                fake_array = []
                for P in range(num_s4):
                    if sp.im(list_s4[P]) ==0:
                        fake_array.append(list_s4[P])
                list_s4 = fake_array
                num_s4 = len(list_s4)

                fake_array2 = list(dict.fromkeys(list_s4))
                list_s4 = fake_array2
                num_s4 = len(list_s4)
            
                fake_array3 = []
                for P in range(num_s4):
                    if  lc_g3.subs(s_4,list_s4[P]) == 0:
                        print("lc(g3)=0!")
                    else:
                        fake_array3.append(list_s4[P])
                list_s4 = fake_array3
                num_s4 = len(list_s4)

                for S in range(num_s4):
                    g3_subs = g3.subs([(s_7,list_s7[Q]),(c_7,list_c7[R]),(s_4,list_s4[S])])
                    deg_g3 = deg_g[3]
                    coef_3 = []
                    for P in range(deg_g3 + 1):
                        n = g3_subs.coeff(c_4, P)
                        n = n.evalf(20)
                        coef_3.insert(0,n)
                    list_c4 = np.roots(coef_3)
                    num_c4 = len(list_c4)
                
                    fake_array = []
                    for P in range(num_c4):
                        if sp.im(list_c4[P]) ==0:
                            fake_array.append(list_c4[P])
                    list_c4 = fake_array
                    num_c4 = len(list_c4)

                    fake_array2 = list(dict.fromkeys(list_c4))
                    list_c4 = fake_array2
                    num_c4 = len(list_c4)
                    
                    for P in range(num_c4):
                        Dict = {}
                        Dict[s_7] = list_s7[Q]
                        Dict[c_7] = list_c7[R]
                        Dict[s_4] = list_s4[S]
                        Dict[c_4] = list_c4[P]
                        Dict[s_1] = 0
                        Dict[c_1] = 1
                        list_Solves.append(Dict)
    else :
        g0 = gb[0].subs([(x,e1),(y,e2),(z,e3)])
        g1 = gb[1].subs([(x,e1),(y,e2),(z,e3)])
        g2 = gb[2].subs([(x,e1),(y,e2),(z,e3)])
        g3 = gb[3].subs([(x,e1),(y,e2),(z,e3)])
        g4 = gb[4].subs([(x,e1),(y,e2),(z,e3)])
        g5 = gb[5].subs([(x,e1),(y,e2),(z,e3)])
        lc_g1 = g1.coeff(c_7,deg_g[1])
        lc_g2 = g2.coeff(s_4,deg_g[2])
        lc_g3 = g3.coeff(c_4,deg_g[3])
        lc_g4 = g4.coeff(s_1,deg_g[4])
        lc_g5 = g5.coeff(c_1,deg_g[5])

        deg_g0 = deg_g[0]
        coef_0 = []

        for W in range(deg_g0+1):
            n = g0.coeff(s_7,W)
            n = n.evalf(20)
            coef_0.insert(0,n)
    
        list_s7 = np.roots(coef_0)
        num_s7 = len(list_s7)
        fake_array = []
        for P in range(num_s7):
            if sp.im(list_s7[P]) ==0:
                fake_array.append(list_s7[P])
        list_s7 = fake_array
        num_s7 = len(list_s7)

        fake_array2 = list(dict.fromkeys(list_s7))
        list_s7 = fake_array2
        num_s7 = len(list_s7)

        fake_array3 = []
        for P in range(num_s7):
            if  lc_g1.subs(s_7,list_s7[P]) == 0:
                print("lc(g1)=0!")
            else:
                fake_array3.append(list_s7[P])
        list_s7 = fake_array3
        num_s7 = len(list_s7)

        
        for Q in range(num_s7):
            g1_subs = g1.subs(s_7, list_s7[Q])
            #list_g1.append(g1_subs)
            deg_g1 = deg_g[1]
            coef_1 = []
            for P in range(deg_g1 + 1):
                n = g1_subs.coeff(c_7, P)
                n = n.evalf(20)
                coef_1.insert(0,n)
            list_c7 = np.roots(coef_1)
            num_c7 = len(list_c7)

            fake_array = []
            for P in range(num_c7):
                if sp.im(list_c7[P]) ==0:
                    fake_array.append(list_c7[P])
            list_c7 = fake_array
            num_c7 = len(list_c7)
        
            fake_array2 = list(dict.fromkeys(list_c7))
            list_c7 = fake_array2
            num_c7 = len(list_c7)
            
            fake_array3 = []
        
            for P in range(num_c7):
                if  lc_g2.subs(c_7,list_c7[P]) == 0:
                    print("lc(g2)=0!")
                else:
                    fake_array3.append(list_c7[P])
        
            list_c7 = fake_array3
            num_c7 = len(list_c7)

            for R in range(num_c7):
                g2_subs = g2.subs([(s_7,list_s7[Q]),(c_7,list_c7[R])])
                deg_g2 = deg_g[2]
                coef_2 = []
                for P in range(deg_g2 + 1):
                    n = g2_subs.coeff(s_4, P)
                    n = n.evalf(20)
                    coef_2.insert(0,n)
                list_s4 = np.roots(coef_2)
                num_s4 = len(list_s4)

                fake_array = []
                for P in range(num_s4):
                    if sp.im(list_s4[P]) ==0:
                        fake_array.append(list_s4[P])
                list_s4 = fake_array
                num_s4 = len(list_s4)

                fake_array2 = list(dict.fromkeys(list_s4))
                list_s4 = fake_array2
                num_s4 = len(list_s4)
                
                fake_array3 = []
                for P in range(num_s4):
                    if  lc_g3.subs(s_4,list_s4[P]) == 0:
                        print("lc(g3)=0!")
                    else:
                        fake_array3.append(list_s4[P])
                list_s4 = fake_array3
                num_s4 = len(list_s4)

                for S in range(num_s4):
                    g3_subs = g3.subs([(s_7,list_s7[Q]),(c_7,list_c7[R]),(s_4,list_s4[S])])
                    deg_g3 = deg_g[3]
                    coef_3 = []
                    for P in range(deg_g3 + 1):
                        n = g3_subs.coeff(c_4, P)
                        n = n.evalf(20)
                        coef_3.insert(0,n)
                    list_c4 = np.roots(coef_3)
                    num_c4 = len(list_c4)
                
                    fake_array = []
                    for P in range(num_c4):
                        if sp.im(list_c4[P]) ==0:
                            fake_array.append(list_c4[P])
                    list_c4 = fake_array
                    num_c4 = len(list_c4)

                    fake_array2 = list(dict.fromkeys(list_c4))
                    list_c4 = fake_array2
                    num_c4 = len(list_c4)

                    fake_array3 = []
                    for P in range(num_c4):
                        if  lc_g4.subs(c_4,list_c4[P]) == 0:
                            print("lc(g4)=0!")
                        else:
                            fake_array3.append(list_c4[P])
                    list_c4 = fake_array3
                    num_c4 = len(list_c4)

                    for T in range(num_c4):
                        g4_subs = g4.subs([(s_7,list_s7[Q]),(c_7,list_c7[R]),(s_4,list_s4[S]),(c_4,list_c4[T])])
                        deg_g4 = deg_g[4]
                        coef_4 = []

                        for P in range(deg_g4 + 1):
                            n = g4_subs.coeff(s_1, P)
                            n = n.evalf(20)
                            coef_4.insert(0,n)
                        list_s1 = np.roots(coef_4)
                        num_s1 = len(list_s1)

                        fake_array = []
                        for P in range(num_s1):
                            if sp.im(list_s1[P]) ==0:
                                fake_array.append(list_s1[P])
                        list_s1 = fake_array
                        num_s1 = len(list_s1)

                        fake_array2 = list(dict.fromkeys(list_s1))
                        list_s1 = fake_array2
                        num_s1 = len(list_s1)

                        fake_array3 = []
                        for P in range(num_s1):
                            if  lc_g5.subs(s_1,list_s1[P]) == 0:
                                print("lc(g5)=0!")
                            else:
                                fake_array3.append(list_s1[P])
                        list_s1 = fake_array3
                        num_s1 = len(list_s1)

                        for U in range(num_s1):
                            g5_subs = g5.subs([(s_7,list_s7[Q]),(c_7,list_c7[R]),(s_4,list_s4[S]),(c_4,list_c4[T]),(s_1,list_s1[U])])
                            deg_g5 = deg_g[5]
                            coef_5 = []

                            for P in range(deg_g5 + 1):
                                n = g5_subs.coeff(c_1, P)
                                n = n.evalf(20)
                                coef_5.insert(0,n)
                            list_c1 = np.roots(coef_5)
                            num_c1 = len(list_c1)

                            fake_array = []
                            for P in range(num_c1):
                                if sp.im(list_c1[P]) ==0:
                                    fake_array.append(list_c1[P])
                            list_c1 = fake_array
                            num_c1 = len(list_c1)

                            fake_array2 = list(dict.fromkeys(list_c1))
                            list_c1 = fake_array2
                            num_c1 = len(list_c1)
                            for P in range(num_c4):
                                Dict = {}
                                Dict[s_7] = list_s7[Q]
                                Dict[c_7] = list_c7[R]
                                Dict[s_4] = list_s4[S]
                                Dict[c_4] = list_c4[T]
                                Dict[s_1] = list_s1[U]
                                Dict[c_1] = list_c1[P]
                                list_Solves.append(Dict)

    TimeSolving = time.time()
    print ("Time for solving the system of algebraic equations:", TimeSolving - TimeDescartes)
    list_SolveTime.append(TimeSolving - TimeDescartes)
    # End of computing roots
    
    # Estimation of errors
    # Get the angle of the joints from s_i and c_i using atan2
    # Solve the forward kinematic problem and calculate the error of the position of the end-effector
    for P in list_Solves :
        S7 = P[s_7]
        C7 = P[c_7]
        S4 = P[s_4]
        C4 = P[c_4]
        S1 = P[s_1]
        C1 = P[c_1]
        theta7 = sp.atan2(S7,C7)
        theta4 = sp.atan2(S4,C4)
        theta1 = sp.atan2(S1,C1)
        t7 = theta7.evalf(10)
        t4 = theta4.evalf(10)
        t1 = theta1.evalf(10)
        if sp.im(theta7.evalf(3)) ==0 and sp.im(theta4.evalf(3)) ==0 and sp.im(theta1.evalf(3)) ==0 :
            print("theta1 =", theta1)
            print("theta4 =", theta4)
            print("theta7 =", theta7)
            X_error = e1 - (-112*cos(t1)*cos(t4)*sin(t7)+16*cos(t1)*cos(t4)-112*cos(t1)*sin(t4)*cos(t7)-136*cos(t1)*sin(t4)+44*sqrt(2)*cos(t1))
            Y_error = e2 - (-112*sin(t1)*cos(t4)*sin(t7)+16*sin(t1)*cos(t4)-112*sin(t1)*sin(t4)*cos(t7)-136*sin(t1)*sin(t4)+44*sqrt(2)*sin(t1))
            Z_error = e3 - (112*cos(t4)*cos(t7)+136*cos(t4)-112*sin(t4)*sin(t7)+16*sin(t4)+104+44*sqrt(2))
            Error = sqrt(X_error**2+Y_error**2+Z_error**2)
            list_Error.append(Error) 
        else :
            print("theta is complex!")
    TimeErrorEstimation = time.time()
    list_ErrorTime.append(TimeErrorEstimation - TimeSolving)
    list_IKSolvingTime.append(TimeErrorEstimation - TimeStart)
    print ("Time for error estimation:", TimeErrorEstimation - TimeSolving)
    print ("Time for solving inverse kinematic problem:", TimeErrorEstimation - TimeStart)

# Shutdown the ox server

libox.ox_shutdown(sv5)
TimeShutdown = time.time()
print ("Shutdown time:", TimeShutdown - TimeErrorEstimation)

# Print the average computing time

print ("-----")

print("Range of sample points: [", min(r), ",", max(r), "]")
print("The average of computing time for choosing the partition:", sum(list_CellNumTime)/len(list_CellNumTime))
print("The average of cpmputing time for the CGS-QE", sum(list_CountRealRootsTime)/len(list_CountRealRootsTime))
print("The average of computing time for substitution", sum(list_SubsTime)/len(list_SubsTime))
print("The average of computing time for the Descartes' rule of signs:", sum(list_DescarteTime)/len(list_DescarteTime))
print("The average of computing time for solving algebraic equations:", sum(list_SolveTime)/len(list_SolveTime))
print("The average of time for inverse kinematics computation:", sum(list_IKSolvingTime)/len(list_IKSolvingTime))
print("The average of computing time for errors:", sum(list_ErrorTime)/len(list_ErrorTime))
print("The average of errors:", sum(list_Error)/len(list_Error))
print("-----")
print("Partition indices used for CGS-QE:", list_CellNum)
print("The number of varid sample points:", len(list_ValidSamplePoint))
print("Valid sample points:", list_ValidSamplePoint)
print("The number of real roots:", list_CountRealRoots)

print (datetime.datetime.now())

exit(0)
