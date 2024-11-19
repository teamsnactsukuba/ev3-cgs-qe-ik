# Inverse kinematic computation with the present method (Rational model)

This section handles the computation of inverse kinematic computation with the present method.

## Computing instructions

### Verification of the shape form

Verification of the shape form can be executed as follows:
```
% asir
load("shape-form-test.rr")$
```
A log is stored in shape-form-test.log.

### Generating the CGS and the partitions data for Python

Generate the CGS and the partitions using ```generate-ab.rr```.
The CGS is stored in ```B.dat``` and the partitions are stored in ```A.dat```.

1. For saving the output, set ``OUTPUT`` variable different from 0, then the output will be saved in ```C.dat```; otherwise, the output will not be saved.
    ```
    % asir
    [2077] OUTPUT = 1$    
    ```
1. Load ```generate-ab.rr``` to execute the computation (the session continues from the above).
    ```
    [2078] load("generate-ab.rr")$
    0
    ```
1. Save the text data of the lists ```AL``` (the partitions) and ```BL``` (the CGS) (defined in ```generate-ab.rr```) as ```A.rr``` and ```B.rr```, respectively.

### Generating the list of degrees of the polynomials in the CGS

Generate the list of degrees of the polynomials in the CGS using ```generate-e.rr```.
The list of the degrees is stored in ```E.dat```.

1. For saving the output, set ``OUTPUT`` variable different from 0, then the output will be saved in ```C.dat```; otherwise, the output will not be saved.
    ```
    % asir
    [2077] OUTPUT = 1$    
    ```
1. Load ```generate-e.rr``` to execute the computation (the session continues from the above).
    ```
    [2078] load("generate-e.rr")$
    0
    ```
1. Save the text data of the lists ```E``` as ```E.rr```.

### Generating sample points for the position of the end-effector 

Sample points for the position of the end-effector have been generated using GenerateSamplePoint.py.

Output is stored in SamplePoint.py, which contains 1000 sample points.

### Inverse kinematic computations

The main program is cgs-qe-ik-N.py (N=0,...,9,10), for the following sample points:
- cgs-qe-ik-0.py: range(0,10) for test
- cgs-qe-ik-N.py, N=1,...,10: range((N-1) * 100, N * 100)

The programs are executed as follows:
```
./cgs-qe-ik-0.py
```
and so on.

Logs are stored in [log](./log/): cgs-qe-ik-N.log for the log of cgs-qe-ik-N.py.

version.log: Version of software used.
