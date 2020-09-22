#!/bin/bash
#gfortran -c fnc.f90
#gfortran -c nr.f90 -L/Users/jamesthorpe/LAPACK/lapack-3.7.0 -llbas -llapack fnc.o
#gfortran extrap.f90 -L/Users/jamesthorpe/LAPACK/lapack-3.7.0 -llbas -llapack
gfortran extrap.f90 -llapack
