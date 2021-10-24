# Inverse kinematic computation with the present method

This section handles the computation of inverse kinematic computation with the present method.

## Computing instructions

### Verification of the shape form

Verification of the shape form can be executed as follows:
```
% asir
load("shape-form-test.rr")$
```
A log is stored in shape-form-test.log.

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
