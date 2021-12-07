include ../Makefile.inc



# --------------------------  SHOULD NOT CHANGE  ------------------------------------------

MKLFLAGS = -I$(MKL_INC) 

FFLAGS = -Ofast -cpp -fopenmp -fPIC

FFLAGS2 = -O2 -cpp -fopenmp -fPIC

CFLAGS = -Ofast  -fPIC

#
# linking to MKL, see /opt/intel/oneapi/mkl/???/documentation/en/common/mkl_link_line_advisor.htm utility
#
# ========== multithreading for mkl =========
#LDFLAGS =  ${F95ROOT}/lib/intel64/libmkl_blas95_lp64.a ${F95ROOT}/lib/intel64/libmkl_lapack95_lp64.a -L${MKLLIB} -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
# ========== sequential for mkl =============
LDFLAGS = ${F95ROOT}/lib/intel64/libmkl_lapack95_lp64.a -L${MKLLIB} -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
# ========== dynamic linking ================
LDFLAGS =    -L${MKLLIB} -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
#LDFLAGS =    -L${MKLLIB} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lgomp -lpthread -lm -ldl
#LDFLAGS =    -L${MKLLIB} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lgomp -lpthread -lm -ldl

# -- F2PY --
F2PYFLAGS = --opt="-Ofast -march=native" --f90flags="-fopenmp -cpp"
export LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64

# ================================================= Nothing to be changed below ==========================================

clean:
	/bin/rm *.o *.mod

.SUFFIXES: .o .f90
.SUFFIXES: .o .f

.c.o:
	$(CC) $(CFLAGS) -c $<

.f90.o: $(INC)
	$(FC) $(FFLAGS) -c $<
# ======================================= librs.a =============================================
mkl_dfti.mod: mkl_dfti.f90
	$(FC) -c $(MKL_INC)/mkl_dfti.f90


# ======================================= PYTHON STUFF ==============================

#for all targets, the most simple. For linux, a1xcorr_python_fc is an alternate way
python: a1xcorr_python


a1xcorPy.o: a1xcorPy.f90
	$(FC) $(FFLAGS) -c -fPIC a1xcorPy.f90

fortranobject.o: 
	$(CC) -I$(PYTHON_INC) -I$(NUMPY_INC) -I$(F2PY_INC) -c -fPIC $(F2PY_SRC) -o fortranobject.o

a1xcorPymodule.c:
	$(F2PY) -m a1xcorPy a1xcorPy.f90

a1xcorPymodule.o: a1xcorPymodule.c
	$(CC)  -I$(PYTHON_INC) -I$(NUMPY_INC) -I$(F2PY_INC) -c -fPIC a1xcorPymodule.c

xcorr.o: xcorr.f90
	$(FC) -c $(MKLFLAGS) $(FFLAGS) xcorr.f90	

a1xcorr_python: mkl_dfti.mod cxcor.o pxcor.o xcorr.o 
	$(F2PY) -c -m _a1xcorPy -DF2PY_REPORT_ON_ARRAY_COPY=1 a1xcorPy.f90 cxcor.o pxcor.o xcorr.o $(LDFLAGS_STATIC); mv _a1xcorPy.*.so ../a1das

a1xcorr_python_fc: mkl_dfti.mod cxcor.o pxcor.o xcorr.o  fortranobject.o a1xcorPy.o a1xcorPymodule.o
	$(FC) -Wall -g -fPIC $(PYTHONFLAG) fortranobject.o a1xcorPymodule.o $(LDFLAGS_STATIC_FC) -lgfortran -o ../a1das/_a1xcorPy.$(SUFFIX)

# ======================================= TEST STUFF (DESTIN) ==============================

test_whitening.o : mkl_dfti.mod  test_whitening.f90
	$(FC) $(FFLAGS) $(MKLFLAGS) -c test_whitening.f90

test_whitening: test_whitening.o 
	$(FC)  $(FFLAGS) $(LDFLAGS) -o $@ mkl_dfti.o  test_whitening.o 

test_bandpass.o : mkl_dfti.mod  test_bandpass.f90
	$(FC) $(FFLAGS) $(MKLFLAGS) -c test_bandpass.f90	

test_bandpass: test_bandpass.o 
	$(FC)  $(FFLAGS) $(LDFLAGS) -o $@ mkl_dfti.o  test_bandpass.o 

test_lowpass.o : mkl_dfti.mod  test_lowpass.f90
	$(FC) $(FFLAGS) $(MKLFLAGS) -c test_lowpass.f90	

test_lowpass: test_lowpass.o 
	$(FC)  $(FFLAGS) $(LDFLAGS) -o $@ mkl_dfti.o  test_lowpass.o 

test: test_whitening test_bandpass test_lowpass 
