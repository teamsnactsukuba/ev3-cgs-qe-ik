# Computing multivariate Hermite quadratic form

This section handles the computation of multivariate Hermite quadratic form and the characteristic polynomial.

1. Note: programs in this directory use ```inconsistent.rr``` in the CGS program by Prof. Katsusuke Nabeshima.
1. For the CGS computed in the previous step, compute the characteristic polynomial of the Hermite quadratic form using Risa/Asir with hermite-compute.rr. Save the output in hermite-compte.log, as follows (the last line ```output()``` should be executed after computation of ``hermite-compute.rr`` is finished):
    ```
    % asir
    output("hermite-compute.log");
    load("hermite-compute.rr");
    output();
    ```
1. Reduce the expression of the formulas in hermite-compute.log with Mathematica. Copy hermite-compute.log to hermite-translate-in.m and 
arrange the contents. Then, use Mathematica to get reduced formula. The result is stored in hermite-translate-out.m.
    ```
    [flopsy5:ev3-cgs-qe-ik/preprocessing-steps/hermite] terui% /Applications/WolframScript.app/Contents/MacOS/wolframscript 
    Wolfram Language 12.0.0 Engine for Mac OS X x86 (64-bit)
    Copyright 1988-2019 Wolfram Research, Inc.  
    In[1]:= << hermite-translate-in.m
    In[2]:= << hermite-translate-out.m
    ```
1. To read the result in Risa/Asir as a list, copy hermite-translate-out.m to hermite-translate-out.rr and change parentheses from {} to [].
1. On Risa/Asir, load hermite-translate-out.rr. The list of characteristic polynomials is defined as C. Save binary data of C for fast loading.
    ```
    % asir
    load("hermite-translate-out.rr")$
    bsave(C, "C.dat")$
    ```
1. In the main steps, the data will be used as
    ```
    % asir
    C = bload("C.dat")$
    ```