# ev3-cgs-qe-ik
An inverse kinematics solver based on the CGS-QE algorithm for an EV3 manipulator

This repository contains an implementation of an inverse kinematics solver based on the CGS-QE algorithm for an EV3 manipulator, described in the following paper:

Shuto Otaki, Akira Terui and Masahiko Mikawa. A Design and an implementation of an inverse kinematics computation in robotics using real quantifier elimination based on Comprehensive Groebner Systems.

## Contents

The contents of this repository are as follows:
- [preprocessing-steps](./preprocessing-steps/): implementation and results of the preprocessing steps.
- [main-steps](./main-steps): implementation and results of the main steps.

## Prerequisites

For executing the implementation here, you need the following:

- Python: version 3 or later
  - NumPy
  - SymPy
- OpenXM infrastructure for communicating mathematical software systems: http://www.openxm.org/
  - Risa/Asir, a computer algebra system: http://www.math.kobe-u.ac.jp/Asir/ (included in the OpenXM distribution and installed automatically with OpenXM under default settings)
- CGS: a program for computing comprehensive Groebner systems in a polynomial ring by Prof. Katsusuke Nabeshima: https://www.rs.tus.ac.jp/~nabeshima/softwares.html
- X Window System server (for executing OpenXM)

Our test environment is as follows:
- Linux 4.15.0 (Ubuntu 18.04)
- Python 3.6.9
- NumPy 1.19.5
- SymPy 1.8
- OpenXM 1.3.3
- Risa/Asir 20210326 (Kobe Distribution)

For executing our implementation, please check the following:
- An environment variable ```OpenXM_HOME``` is set properly
- The following library files exist in ```${OpenXM_Home}/lib```:
  - For UNIX/Linux systems:
    - libgmp.so
    - libgc.so
    - libmpfr.so
    - libox.so
  - For macos:
    - libgmp.dylib
    - libgc.dylib
    - libmpfr.dylib
    - libox.dylib
- Python's ```ctypes``` package works properly