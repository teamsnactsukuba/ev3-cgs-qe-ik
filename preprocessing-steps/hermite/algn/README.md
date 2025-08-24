# Computing multivariate Hermite quadratic form (algebraic number version)

This section handles the computation of multivariate Hermite quadratic form and the characteristic polynomial.

For the CGS computed in the previous step, the programs in this directory is used to compute the characteristic polynomial of the Hermite quadratic form using Risa/Asir.

Note: programs in this directory use ```inconsistent.rr``` in the CGS program by Prof. Katsusuke Nabeshima. See: https://www.rs.tus.ac.jp/~nabeshima/softwares.html

## Generating the Hermite quadratic form for each segment

 For saving the output, set ``OUTPUT`` variable different from 0, then the output will be saved in ```C.dat```; otherwise, the output will not be saved.
    ```
    % asir
    OUTPUT = 1$    
    ```
1. Load ```hermote-compute.rr```.  Then, for the segment (S_i, G_i), where S_i is the defining polynomials for the segment and G_i is its Groebner basis, it computes the quadratic form of the residue class ring Q[x_1,...,x_n]/<G_i>.
execute the computation (the session continues from the above).
If ```OUTPUT = 1```, the quadratic form of each S_i is saved as "H-i.dat".
    ```
    load("hermite-compute.rr")$
    0
    ```
1. Or, if you want to generate the quadratic form of the residue class ring separately, you can do as follows:
    ```
    % asir
    load("hermite-compute-i.rr")$
    ```
1. Then, if you want to generate the quadratic form of the residue class ring corresponding to S_i, 
    ```
    hermite_i(I, Mode)
    ```
    * For ```Mode = 1```, it calculates the quadratic form w.r.t. variables c_1,s_1,c_4,s_4,c_7,s_7. 
    On the other hand, for ```Mode = 2```, it calculates the quadratic form w.r.t. variables c_4,s_4,c_7,s_7 (c_1 and s_1 are treated as free variables). 
    * If ```OUTPUT = 1```, the quadratic form of each S_i is saved as "H-i.dat".

## Arranging the quadratic forms to pack into one set

For this purpose, use ```hermite-generate-c.rr```.

    ```
    % asir
    load("hermite-generate-c.rr")$
    ```
This program does the following.

1. It reads ```H-i.dat```.
1. Multiply the denominator to make ```Hi``` a polynomial.
1. Calculate the remainder of Hp divided by a^2-2 which corresponds a = 2^(1/2).
1. Substitute ```a``` with ```2^(1/2).```
1. Collect the quadratic forms as ```C``` and, if ```OUTPUT = 1```, save ```C``` to　```C.dat```.

```C.dat``` is used in the main step.



