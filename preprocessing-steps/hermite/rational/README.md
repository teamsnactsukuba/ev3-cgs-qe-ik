# Computing multivariate Hermite quadratic form (Rational version)

This section handles the computation of multivariate Hermite quadratic form and the characteristic polynomial.

1. Note: programs in this directory use ```inconsistent.rr``` in the CGS program by Prof. Katsusuke Nabeshima. See: https://www.rs.tus.ac.jp/~nabeshima/softwares.html
1. For the CGS computed in the previous step, compute the characteristic polynomial of the Hermite quadratic form using Risa/Asir with hermite-compute.rr. For saving the output, set ``OUTPUT`` variable different from 0, then the output will be saved in ```C.dat```; otherwise, the output will not be saved.
    ```
    % asir
    [2077] OUTPUT = 1$    
    ```
1. Load ```hermote-compute.rr``` to execute the computation (the session continues from the above).
    ```
    [2078] load("hermite-compute.rr")$
    0
    ```
1. To see the contents, load the data as follows:
    ```
    % asir
    C = bload("C.dat")$
    ```