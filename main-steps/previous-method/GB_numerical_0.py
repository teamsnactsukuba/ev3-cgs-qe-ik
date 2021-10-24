#!/usr/bin/env python3
# coding: UTF-8   
# An inverse kimatics solver for EV3 with solving univariate equations
# Akira Terui, December 21, 2019: the original program
# Noriyuki Horigome, January 2020: for experiment of his master's thesis
# Shuto Otaki, December 30, 2020: with new EV3 model
# Masahiko Mikawa, July 2021: ported to Python 3
# Akira Terui, August 2021: code refactoring

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
import sys
import platform
import os
import SamplePoint
import ev3_model

SamplingPoint = []
GBTimeList = []
PrintTimeList = []
SolveTimeList = []
ShutdownTimeList = []
IKTimeList = []
ErrorTimeList = []
ErrorList = []
RelativeErrorList = []

t1 = ev3_model.t1
t4 = ev3_model.t4
t7 = ev3_model.t7
x = Symbol('x')
s1 = Symbol('s1')
s4 = Symbol('s4')
s7 = Symbol('s7')
c1 = Symbol('c1')
c4 = Symbol('c4')
c7 = Symbol('c7')

M = ev3_model.M

e1 = M[0,3]
e2 = M[1,3]
e3 = M[2,3]

# Calculate substitution of trigonometric functions by variables

trigSubs = [
    (sin(t1), s1),
    (sin(t4), s4), 
    (sin(t7), s7),
    (cos(t1), c1),
    (cos(t4), c4),
    (cos(t7), c7)
    ]

e11 = e1.subs(trigSubs)
e21 = e2.subs(trigSubs)
e31 = e3.subs(trigSubs)

# print("e1 =", e1)
# print("e2 =", e2)
# print("e3 =", e3)

# print()

# print("e11 =", e11)
# print("e21 =", e21)
# print("e31 =", e31)

# print()

#--------------------------------------------------------------------------------------------------------------------------------------------------------

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

# Start the ox server

sv5 = libox.ox_start(OX_SERVER_HOST, OX_PROG1, OX_PROG2)

r = range(0,10) # range of sample points

for k in r:

    # Take a point from SamplePoint.py
    # Let x,y,z be each coordinates
    Points0 = SamplePoint.Points[k]
    Point = Points0
    X_Coordinate = Point[0]
    Y_Coordinate = Point[1]
    Z_Coordinate = Point[2]
    
    Time1 = time.time()
    
    # Construct generators of the ideal

    f = []

    f.append(str(e11) + " - (" + str(X_Coordinate) + ")")
    f.append(str(e21) + " - (" + str(Y_Coordinate) + ")")
    f.append(str(e31) + " - (" + str(Z_Coordinate) + ")")
    f.append("c1*c1 + s1*s1 - 1")
    f.append("c4*c4 + s4*s4 - 1")
    f.append("c7*c7 + s7*s7 - 1")
    term_order_string = "[c1,s1,c4,s4,c7,s7]"
    for i in range(len(f)):
        f[i] = f[i].replace("sqrt(2)", "2^(1/2)")
        #print "f[" + str(i) + "] = " + f[i]
    
    print("Case:", k)
    print("term_order:", term_order_string)
    
    # Construct command string for Asir

    ox_command_string = "G = tolex_tl(gr([" + f[0] + ", " + f[1] + ", " + f[2] + ", " + f[3] + ", " + f[4] + ", " + f[5] + "], " + term_order_string + ", 0),"+str(term_order_string)+",2,"+str(term_order_string)+",1)$"
    print("ox_command_string:", ox_command_string)
    
    OX_COMMAND_GB = create_string_buffer(ox_command_string.encode('utf-8'))

    # Compute Groebner basis with Asir via OpenXM

    libox.ox_execute_string(sv5, OX_COMMAND_GB)

    # time.sleep(1.5)

    # Extract the size of G-base
    # for G-base = [g_1, ... , g_m] , let gbase_length := m

    ox_command_string = "length(G);"
    OX_COMMAND = create_string_buffer(ox_command_string.encode('utf-8'))
    libox.ox_execute_string(sv5, OX_COMMAND)
    OX_OUTPUT = libox.ox_popString(sv5)
    gbase_length = int(OX_OUTPUT)

    # Extract each polynomial in G-base one-by-one

    GB = []

    for i in range(gbase_length):
        ox_command_string = "G[" + str(i) + "];"
        OX_COMMAND = create_string_buffer(ox_command_string.encode('utf-8'))
        libox.ox_execute_string(sv5, OX_COMMAND)
        OX_OUTPUT = libox.ox_popString(sv5)
        f = OX_OUTPUT.decode("utf-8")
        f = f.replace('^','**')
        GB.append(sympify(f))

    Time2 = time.time()
    GBTime = Time2 - Time1
    print("GBTime:", GBTime)
    GBTimeList.append(GBTime)

    print("GB =", GB)
    point = [X_Coordinate, Y_Coordinate, Z_Coordinate]
    SamplingPoint.append(k)
    print("X =", X_Coordinate)
    print("Y =", Y_Coordinate)
    print("Z =", Z_Coordinate)
    
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Solving univariate equations

    Time3 = time.time()
    PrintTime = Time3 - Time2
    print("PrintTime:", PrintTime)
    PrintTimeList.append(PrintTime)

    g0 = GB[0]
    g1 = GB[1]
    g2 = GB[2]
    g3 = GB[3]
    g4 = GB[4]
    g5 = GB[5]
    g1_list = []

    # The leading terms of g1,...,g5
    LCP_g1 = g1.coeff(c7, degree(g1, gen=c7))
    LCP_g2 = g2.coeff(s4, degree(g2, gen=s4))
    LCP_g3 = g3.coeff(c4, degree(g3, gen=c4))
    LCP_g4 = g4.coeff(s1, degree(g4, gen=s1))
    LCP_g5 = g5.coeff(c1, degree(g5, gen=c1))
        
    deg_poly = polys.degree(g0)
    coefficient_1 = []
    for i in range(deg_poly+1):
        a = g0.coeff(s7,i)
        a = a.evalf(20)
        coefficient_1.insert(0,a)
    # Solve g0 == 0 for s7
    s7_list = np.roots(coefficient_1)

    num_s7 = len(s7_list)

    fake_array = []
    for p in range(num_s7):
        if im(s7_list[p]) == 0:
            fake_array.append(s7_list[p])

    s7_list = fake_array
    num_s7 = len(s7_list)

    # Test if the partial solutions satisfy the Extension Theorem
    fake_array2 = []
    for p in range(num_s7):
        if  LCP_g1.subs(s7,s7_list[p]) == 0:
            print("Makes the leading term 0:", s7_list[p])
        else:
            fake_array2.append(s7_list[p])
    s7_list = fake_array2
    num_s7 = len(s7_list)

    Solves = []
    for i in range(num_s7):  # Solve for the other solutions for s7
        g1subs = g1.subs(s7, s7_list[i])
        g1_list.append(g1subs)
        
        deg_poly = polys.degree(g1subs)
        coefficient_2 = []

        for j in range(deg_poly+1):
            b = g1subs.coeff(c7, j)
            b = b.evalf(20)
            coefficient_2.insert(0,b)

        # Solve g1 == 0 for c7
        c7_list = np.roots(coefficient_2)

        num_c7 = len(c7_list)

        # Test if the partial solutions satisfy the Extension Theorem
        fake_array2 = []
        for p in range(num_c7):
            if  LCP_g2.subs(c7, c7_list[p]) == 0:
                print("Makes the leading term 0:", c7_list[p])
            else:
                fake_array2.append(c7_list[p])
        c7_list = fake_array2
        num_c7 = len(c7_list)
        
        for j in range(num_c7): # Solve for the other solutions for s7 and c7
            g2subs = g2.subs([(s7,s7_list[i]),(c7,c7_list[j])])
            deg_poly = polys.degree(g2subs)

            coefficient_3 = []
            for k in range(deg_poly + 1):
                c = g2subs.coeff(s4,k)
                c = c.evalf(20)
                coefficient_3.insert(0,c)

            # Solve g2 == 0 for s4            
            s4_list = np.roots(coefficient_3)
            
            num_s4 = len(s4_list)

            # Test if the partial solutions satisfy the Extension Theorem
            fake_array2 = []
            for p in range(num_s4):
                if  LCP_g3.subs(s4, s4_list[p]) == 0:
                    print("Makes the leading term 0:", s4_list[p])
                else:
                    fake_array2.append(s4_list[p])
            s4_list = fake_array2
            num_s4 = len(s4_list)


            for k in range(num_s4): # Solve for the other solutions for s7, c7 and s4
                g3subs = g3.subs([(s4,s4_list[k]),(c7,c7_list[j]),(s7,s7_list[i])])
                deg_poly = polys.degree(g3subs)

                coefficient_4 = []
                for l in range(deg_poly + 1):
                    d = g3subs.coeff(c4,l)
                    d = d.evalf(20)
                    coefficient_4.insert(0,d)

                # Solve g3 == 0 for c4
                c4_list = np.roots(coefficient_4)

                num_c4 = len(c4_list)

                # Test if the partial solutions satisfy the Extension Theorem
                fake_array2 = []
                for p in range(num_c4):
                    if  LCP_g4.subs(c4,c4_list[p]) == 0:
                        print("Makes the leading term 0:", c4_list[p])
                    else:
                        fake_array2.append(c4_list[p])
                c2_list = fake_array2
                num_c5 = len(c4_list)


                for l in range(num_c4): # Solve for the other solutions for s7, c7, s4 and c4
                    g4subs = g4.subs([(c4,c4_list[l]),(s4,s4_list[k]),(c7,c7_list[j]),(s7,s7_list[i])])
                    deg_poly = polys.degree(g4subs)
                    coefficient_5 = []

                    for m in range(deg_poly + 1):
                        e = g4subs.coeff(s1,m)
                        e = e.evalf(20)
                        coefficient_5.insert(0,e)

                    # Solve g4 == 0 for s1
                    s1_list = np.roots(coefficient_5)

                    num_s1 = len(s1_list)

                    # Test if the partial solutions satisfy the Extension Theorem
                    fake_array2 = []
                    for p in range(num_s1):
                        if  LCP_g5.subs(s1,s1_list[p])==0:
                            print("Makes the leading term 0:")
                        else:
                            fake_array2.append(s1_list[p])
                    s1_list = fake_array2
                    num_s1 = len(s1_list)
                    
                    
                    for m in range(num_s1): # Solve for the other solutions for s7, c7, s4, c4 and s1
                        g5subs = g5.subs([(s1,s1_list[m]),(c4,c4_list[l]),(s4,s4_list[k]),(c7,c7_list[j]),(s7,s7_list[i])])
                        deg_poly = polys.degree(g5subs)

                        coefficient_6 = []
                        for n in range(deg_poly +1):
                            f = g5subs.coeff(c1,n)
                            f = f.evalf(20)
                            coefficient_6.insert(0,f)

                        # Solve g5 == 0 for c1                        
                        c1_list = np.roots(coefficient_6)
                        num_c1 = len(c1_list)
                        
                        for n in range(num_c1):
                            Dict = {}
                            Dict[s7] = s7_list[i]
                            Dict[c7] = c7_list[j]
                            Dict[s4] = s4_list[k]
                            Dict[c4] = c4_list[l]
                            Dict[s1] = s1_list[m]
                            Dict[c1] = c1_list[n]
                            Solves.append(Dict)

    print("Soltions:", Solves)

    Time4 = time.time()
    SolvingTime = Time4 - Time3
    print("Solving Time:", SolvingTime)
    SolveTimeList.append(SolvingTime)
    
    IKSolvingTime = Time4 - Time1
    print("IKSolvingTime:", IKSolvingTime)
    IKTimeList.append(IKSolvingTime)
    
    for Trigo in Solves:
        # print(Trigo)
        S1 = Trigo[s1]
        C1 = Trigo[c1]
        S4 = Trigo[s4]
        C4 = Trigo[c4]
        S7 = Trigo[s7]
        C7 = Trigo[c7]

        theta1 = atan2(S1,C1)
        theta4 = atan2(S4,C4)
        theta7 = atan2(S7,C7)

        if im(theta1.evalf(3)) == 0 and im(theta4.evalf(3)) == 0 and im(theta7.evalf(3)) == 0: 
            print("theta1 =", theta1)
            print("theta2 =", theta4)
            print("theta3 =", theta7)
            print("----------------------")
            X_error = e1.subs([(t1,theta1.evalf(10)),(t4,theta4.evalf(10)),(t7,theta7.evalf(10))]) - X_Coordinate
            Y_error = e2.subs([(t1,theta1.evalf(10)),(t4,theta4.evalf(10)),(t7,theta7.evalf(10))]) - Y_Coordinate
            Z_error = e3.subs([(t1,theta1.evalf(10)),(t4,theta4.evalf(10)),(t7,theta7.evalf(10))]) - Z_Coordinate
            # X_error = f1.subs([(t1,(atan2(S1,C1)).evalf(10)),(t2,(atan2(S2,C2)).evalf(10)),(t3,(atan2(S3,C3)).evalf(10))]) - X_Coordinate
            # Y_error = f2.subs([(t1,(atan2(S1,C1)).evalf(10)),(t2,(atan2(S2,C2)).evalf(10)),(t3,(atan2(S3,C3)).evalf(10))]) - Y_Coordinate
            # Z_error = f3.subs([(t1,(atan2(S1,C1)).evalf(10)),(t2,(atan2(S2,C2)).evalf(10)),(t3,(atan2(S3,C3)).evalf(10))]) - Z_Coordinate
            Error = sqrt(X_error**2 + Y_error**2 + Z_error**2)
            ErrorList.append(Error)
            print("----------------------")
        else:
            print("Oops, theta(s) are complex number!")
       
    Time5 = time.time()
    ErrorCalculatingTime = Time5 - Time4
    print("ErrorCaluclatingTime:", ErrorCalculatingTime)
    ErrorTimeList.append(ErrorCalculatingTime)

    # IKTimeList.append(Time5-Time1)

    print("----------------------")

# Shutdown the ox server

libox.ox_shutdown(sv5)

# Print the average computing time

print("Range of sample points: [", min(r), ",", max(r), "]")
print("The average of computing time for Groebner basis (Time2 - Time1):", sum(GBTimeList)/len(GBTimeList))
print("The average of printing time (Time3 - Time2):", sum(PrintTimeList)/len(PrintTimeList))
print("The average of computing time for solving algebraic equiations (Time4 - Time3):", sum(SolveTimeList)/len(SolveTimeList))
print("The average of time for inverse kinematics computation (Time4 - Time1):", sum(IKTimeList)/len(IKTimeList))
print("The average of computing time for errors (Time5 - Time4):", sum(ErrorTimeList)/len(ErrorTimeList))
print("The average of errors:", (sum(ErrorList)/len(ErrorList)).evalf(10))
print("Sample points:", SamplingPoint)
print("The number of sample points:", len(SamplingPoint))

sys.exit(0)



