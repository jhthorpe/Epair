#!/bin/bash
gfortran -c jac.f90 -llapack /Users/james.thorpe/BLAS-3.8.0/*.a 
gfortran -c fmin.f90 -llapack /Users/james.thorpe/BLAS-3.8.0/*.a
gfortran -c nlslv.f90 -llapack /Users/james.thorpe/BLAS-3.8.0/*.a
gfortran extrap.f90 -llapack /Users/james.thorpe/BLAS-3.8.0/*.a *.o
