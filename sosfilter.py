# sosfilter module
# Wrapper for shared memory version for the sos filtering of an array of signal
# Syntax is identical to sosfilt and sosfiltfilt from scipy except that axis of filtering is
# only along 1st axis.
#
# A fortran module including openMP parallelized version of the filtering must be compiled
# using f2py utility
#
# Under linux, compilation using gfortran is :
#   f2py -c -m sosfilterOMP  -DF2PY_REPORT_ON_ARRAY_COPY=1 -lgomp --opt="-Ofast -march=native" --f90flags="-fopenmp -cpp" sosfilteringOMP.f90
# With intel compiler it becomes:
#   f2py -c -m sosfilterOMP  -DF2PY_REPORT_ON_ARRAY_COPY=1 -L/opt/intel/compilers_and_libraries_2017.4.181/mac/compiler/lib -liomp5 --f90flags="-qopenmp -fpp" sosfilteringOMP.f90
#
# The filtering can be applied to single/double precision array of data

import sosfilterOMP

def sosfilt(sos,sig):
    '''
    sosfilt(sos,sig)
    Apply a sos filter along the fast second dimension of the 2D array sig
    using a openmp parallelization scheme
    sos: (norder, ncof) matrix of sos coefficients
    sig: (nsignal, npoint) 2D array of signal(s) to be filtered along axis 1
    '''
    if (sig.dtype=='float64'):
       sosfilterOMP.sosfilt_d(sos,sig)
    elif (sig.dtype=='float32'):
       sosfilterOMP.sosfilt_s(sos,sig)
    else:
       printf('sosfiltOMP not implemented for ',sig.dtype)
       exit()

    return sig

def sosfiltfilt(sos,sig):
    '''
    sosfiltfilt(sos,sig)
    Apply a sos filter along the fast second dimension of the 2D array sig
    using a openmp parallelization scheme, along direct and backward direction (null phase filter)
    sos: (norder, ncof) matrix of sos coefficients
    sig: (nsignal, npoint) 2D array of signal(s) to be filtered along axis 1
    '''
    if (sig.dtype=='float64'):
       sosfilterOMP.sosfiltfilt_d(sos,sig)
    elif (sig.dtype=='float32'):
       sosfilterOMP.sosfiltfilt_s(sos,sig)
    else:
       printf('sosfiltOMP not implemented for ',sig.dtype)
       exit()

    return sig
