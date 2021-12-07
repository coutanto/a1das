# Some dynamic python object that uses the MKL libraries needs to be compiled with the static libraries
# this requires under linux to use particular link option that cannot be passed as arguments
# For those modules, we nned to create the python binary module without f2py
# 1) We first create  the c interface "raw2strPymodule.c" using f2py -m
# 2) then make the link
# building with numpy > 1.18.2 is broken (1.20, 1.21) due to unresolve symbol npy_PyErr_ChainExceptionsCause()
#
#
BASE = $(PWD)/..

# --------------------------  TO BE SET ACCORDING TO YOUR COMPUTER  ----------
#include ../Makefile.linux
#include ../Makefile.luke
include $(BASE)/Makefile.mac


# --------------------------  TO BE SET ACCORDING TO YOUR COMPUTER (LESS COMMON) ----------
# if changing compiler, change f2py compiler options accordingly: FFLAGS, F2PYFLAGS
FC=gfortran
CC=gcc
F2PY = f2py
