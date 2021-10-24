# Computing Comprehensive Groebner Systems (CGS) for inverse kinematic computation

This section handles the computation of Comprehensive Groebner Systems (CGS).

## Computing instructons

1. Programs in this directory use CGS program by Prof. Katsusuke Nabeshima.
1. Computation of the CGS of the basis F can be executed as follows:
    ```
    % asir
    load("cgs-F-input.rr")$
    ```
1. Computation of the CGS of the basis H can be executed as follows:
    ```
    % asir
    load("cgs-H-input.rr")$
    ```
## Contents

### CGS

The CGS and related data are stored as follows:

- F.rr: the CGS of F (in text format)
- F.dat: the CGS of F (in binary format)
- F-basis contains Groebner bases of each segment of the CGS of F
  - G-NN.rr: the Groebner basis in segment NN
- F-segments contains the defining polynomials of each segment of the CGS of F
  - F-NN-1.rr and F-NN-2.rr: defining polynomials of the segment NN in the form of V(F-NN-1) \ V(F-NN-2), where V denotes the affine variety
- H.rr: the CGS of H (in text format)
- H.dat: the CGS of H (in binary format)
- H-basis contains Groebner bases of each segment of the CGS of H
  - H-NN.rr: the Groebner basis in segment NN
- H-segments contains the defining polynomials of each segment of the CGS of H
  - H-NN-1.rr and H-NN-2.rr: defining polynomials of the segment NN in the form of V(H-NN-1) \ V(H-NN-2), where V denotes the affine variety

### Verification of the existence of real points 

Verification of the existence of real points was executed by hand using Risa/Asir and Mathematica. 

Logs of the verification are stored as follows:

- The segments of the CGS of F: see 
[F-segments/REAEME.md](./F-segments/REAEME.md) 
- The segments of the CGS of H: see 
[H-segments/REAEME.md](./H-segments/REAEME.md) 


