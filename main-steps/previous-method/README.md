# Inverse kinematic computation with the previous method

This section handles the computation of inverse kinematic computation with the previous method describe in:

N. Horigome, A. Terui, M. Mikawa. A Design and an Implementation of an Inverse Kinematics Computation in Robotics Using Gröbner Bases. Proceedings of the 7th International Congress on Mathematical Software (ICMS 2020). Lecture Notes in Computer Science 12097, Springer, 2020, 3–13. https://doi.org/10.1007/978-3-030-52200-1_1

## Computing instructions

### Generating sample points for the position of the end-effector 

Sample points for the position of the end-effector have been generated using GenerateSamplePoint.py in [present method](../present-method).

Output is stored in SamplePoint.py, which contains 1000 sample points.

### Inverse kinematic computations

The main program is GB_numerical_N.py (N=0,...,9,10), for the following sample points:
- GB_numerical_0.py: range(0,10) for test
- GB_numerical_N.py, N=1,...,10: range((N-1) * 100, N * 100)

The programs are executed as follows:
```
./GB_numerical_0.py
```
and so on.

Logs are stored in [log](./log/): GB_numerical_N.log for the log of GB_numerical_N.py.

version.log: Version of software used.

host-environment.log: Specs of the host environment.
