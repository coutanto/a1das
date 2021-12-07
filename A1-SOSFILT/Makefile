include ../Makefile.inc

# --------------------------  SHOULD NOT CHANGE  ------------------------------------------

FFLAGS = -Ofast -cpp -fopenmp -fPIC

FFLAGS2 = -O2 -cpp -fopenmp -fPIC

CFLAGS = -Ofast  -fPIC

#
# -- F2PY --
F2PYFLAGS = --opt="-Ofast -march=native" --f90flags="-fopenmp -cpp" -DF2PY_REPORT_ON_ARRAY_COPY=1
# ================================================= Nothing to be changed below ==========================================

clean:
	/bin/rm *.o *.mod 

.SUFFIXES: .o .f90
.SUFFIXES: .o .f

.c.o:
	$(CC) $(CFLAGS) -c $<

.f90.o: $(INC)
	$(FC) $(FFLAGS) -c $<
# ======================================= PYTHON STUFF ==============================
python: sosfilter

sosfilter:
	f2py -c -m _sosfilterPy $(F2PYFLAGS) -lgomp sosfilteringOMP.f90 ; mv _sosfilterPy.*.so ../a1das


