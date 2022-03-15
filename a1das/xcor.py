#
# Module to compute xcorr from strain(rate) binary data read from files of
# real-time DAS data-flux
#
# O. Coutant, ISTerre, UGA, 2021
# D. Nziengui, Febus, 2021
#
# Content:
#  Interface user - binary (fortran) to compute cross-correlation from a Febus file/socket_stream
#
#    register_couple() = register a set of couple of trace {i,j} for which we want to compute xcorr
#    register_par()    = register the process parameters (onebit, whitening etc)
#    reset_couple()    = reset the list of trace couple {i,j}
#    trace_index_list()= give the list of trace involved in computing xcor
#    get_nxcorr()      = return the number of cross-correlations registered
#    compute_xcorr()   = compute the cross-correlations registered
#
#  use binary module _a1xcorPy.???.so
#
#



import numpy as np

__version__ = "a1xcor Version 1.0.1"
__doc__='Functions to compute cross-correlations. Optimization with fortran codes using openMP and MKL libraries'

def register_couple(couple, ntrace, base=0):
    """
    ## Description
    Register the list of trace couples {i,j} for with we compute the cross-correlation

    if base = 0 (default = python, C convention)
        0 <= i < ntrace
        0 <= j < ntrace

    if base = 1 (Fortran, matlab convention)
        1 <= i <= ntrace
        1 <= j <= ntrace

    ## Input:
        couple= 2D ndarray (or list, or tuple) of size [ncouple , 2] (ex: [[1,1],[1,2],[1,3]] to compute {1,1}, {1,2}, {1,3} cross-correlations
        trace= number of traces in the 2D das section
        base = starting index for trace number, 0 (default) or 1

    ## Return:
        status= 0 ok; 1 wrong index number (<0 or >ntrace)

    ##Usage example
        >>>#For a section containing 500 traces
        >>>#Compute correlation of trace 0 and 1 with all other traces
        >>>import a1das
        >>>ntrace=500
        >>>couple =  [list([0,j]) for j in range(0,ntrace)]
        >>>couple += [list([1,j]) for j in range(0,ntrace)]
        >>>a1das.xcor.register_couple(couple,ntrace=ntrace)


    """

    # import the fortran built _a1xcorPy module
    from a1das import _a1xcorPy
    from ._a1das_exception import DataTypeError, WrongValueError

    if isinstance(couple, list) or isinstance(couple, tuple):
        try:
            couple = np.array(couple)
        except:
            DataTypeError('could not convert couple list into a numpy ndarray of dimension ncouple x 2')
    elif not isinstance(couple, np.ndarray):
        raise DataTypeError('couple list must be a numpy ndarray of dimension ncouple x 2 (or a list or tuple)')
    #if not isinstance(couple, np.ndarray):
    #    raise DataTypeError('couple list must be a numpy ndarray of dimension ncouple x 2')
    if couple.shape[1] != 2:
        raise DataTypeError('couple list must be a numpy ndarray of dimension ncouple x 2')

    if base == 0:
        couple = couple + 1

    status = _a1xcorPy.a1xcor.register_couple(couple.astype('int32'), ntrace)
    if status > 0:
        raise WrongValueError('one set of couple indices is wrong')

    if base == 0:
        couple = couple - 1


    return status

def register_par(proc_key, val=[]):
    """
    Register the parameters used to pre-process data for correlation computation
    proc_key: string referring to a known process

    """
    from a1das import _a1xcorPy
    from ._a1das_exception import DataTypeError, WrongValueError

    if not isinstance(proc_key, str):
        raise DataTypeError('process key must be a string referring to a processing')
    status = _a1xcorPy.a1xcor.register_par(proc_key, np.array(val))
    if status == 1:
        _a1xcorPy.a1xcor.display_pxc()
        raise WrongValueError('Wrong processing name <'+proc_key+'>')
    if status == 2:
        raise WrongValueError('Wrong number of argument for processing <'+proc_key+'>')

def list_processing():
    """
    ##Description
    print the list of pre-processing available when computing cross-correlation
    """
    from a1das import _a1xcorPy
    _a1xcorPy.a1xcor.display_pxc()


def reset_couple():
    """
    Reset the list of couple of trace registered for xcorr computation
    """
    from a1das import _a1xcorPy
    _a1xcorPy.a1xcor.reset_couple()

def trace_index_list():
    """
    trace_index_list: return the indices of the traces used to compute correlation
    """

    from a1das import _a1xcorPy
    #
    # retrieve the number of trace currently declared
    #
    ntrace = _a1xcorPy.a1xcor.get_ntrace()
    if ntrace == 0:
        print('No traces/couple_of_traces registered yet')
        return None
    #
    # read the list of trace used
    #
    list, nlist = _a1xcorPy.a1xcor.trace_index_list(np.zeros((ntrace,), dtype='int32'))

    list = np.resize(list, (nlist,))

    return list

def compute_xcorr(arrayIn, lag=None, stack=False, verbose=0, base=0):
    """
    ## Description
    xcorr, lags, ijx, ier = compute_xcorr(arrayIn, lag=None, stack=False, verbose=0)

    Compute cross-correlations in a 2D array of traces [Ntime x Nspace]
    The array is given either as a 2D ndarray or as an A1Section class instance.

    Cross-correlations are computed for the couple of traces {i,j} declared previously by a1das.xcor.register_couple()
    and using the parameter declared previously by a1das.xcor.register_par()

    !!!!! WARNING !!!! currently process only float64 data array
    
    TODO: write a float32 version, write a transposed version accepting array[nspace x ntime]

    ## Input
        arrayIn: 2D (float64) ndarray [nTime x nSpace], or <A1Section> class instance (with float64 data array)
        lag:    (float) time lag, in sample for ndarray input, in sec for <A1Section>
        stack:  (bool) stack the correlations for consecutive windows of 2*lag+1 length if 2*lag+1 < length signal (default False)
        base:   (int) <0> if first index of trace is 0, <1> if first index is 1 (default 0)
        verbose: verbosity (default 0)

    ## Return
        xcorr:   2D ndarray of cross-correlations [nxcor x nTime] (float64)
        lags:    vector of lag times, in sample for ndarray input, in sec for <A1Section> input
        ijx:     2D ndarray of size[nxcor, 2] giving trace index for each correlation
        ier: return code, 0 if ok

    ## Usage example
        >>> import a1das
        >>> f = a1das.open(filename, format='febus')
        >>> a1 = f.read()
        >>> a1das.xcor.register_couple([[1,2],[1,3],[1,4]]) # compute xcorr between traces 1&2, 1&3, 1&4
        >>> a1das.xcor.register_par('onebit') # set onebit pre-processing
        >>> lag = 5.           # set lag time to 5sec
        >>> # compute xcorr between -5sec and 5sec, compute FFT over consecutive windows of 10 sec and stack them
        >>> xcorr, lags, ijx, ier = compute_xcorr(a1, lag=lag, stack = True)
        >>> # same but for ndarray with no stack,
        >>> xcorr, lags, ijx, ier = compute_xcorr(a1.data, lag=5./a1['dt'], stack = False)
    """
    from a1das import _a1xcorPy, A1Section
    from a1das._a1das_exception import DataTypeError
    from numpy import ndarray

    #
    # determine the cross-correlation length depending on lag and parity of signal length
    #
    #                    xcor_length    lag            mkl_fft_shift        copy_indices
    #  ntime even(2*k)
    #   lag = 0             ntime       ntime/2           ntime/2 (=lag)      [1:ntime]
    #   lag # 0            2*lag+1        --              lag                 [1:2*lag+1]
    # ntime odd (2*k+1)
    #   lag = 0              ntime      (ntime-1)/2      (ntime-1)/2(=lag)    [1:ntime](=[1:2*lag+1])
    #   lag # 0             2*lag+1        --            lag                  [1:2*lag+1]
    #
    #    cross-correlation zero lag is thus always at position <lag> in the xcorr array (xcor[lag-1])
    #

    lag_is_second = False
    if isinstance(arrayIn, A1Section):
        a1 = arrayIn
        arrayIn = arrayIn.data
        lag = int(lag/a1.data_header['dt'])
        if a1.data_header['axis1'] == 'space':
            transposed = True
            ntime = arrayIn.shape[1]
        else:
            transposed = False
            ntime = arrayIn.shape[0]
        lag_is_second = True
        if arrayIn.dtype != np.float64:
            raise DataTypeError('Can only process float64 data type, please convert before calling xcorr')
    elif isinstance(arrayIn, ndarray):
        ntime = arrayIn.shape[0]
        transposed = False
    else:
        raise DataTypeError('Data array must be an A1Section or a numpy ndarray (ntime x nspace)')

    if lag is None:
        lag = ntime
    lag = int(lag)
    if ntime % 2 == 0: # Even case ntime=2*k
        if lag is None or (lag == 0) or (2*lag + 1 > ntime):
            lag = ntime/2
            len = ntime
            lags = np.arange(-lag,lag)
        else:
            len = 2*lag+1
            lags = np.arange(-lag,lag+1)

    else:             # Odd case ntime=2*k+1
        if lag is None or lag == 0:
            lag = (ntime-1)/2
            len = ntime
        else:
            if (2*lag + 1 > ntime):
                lag = (ntime-1)/2
            len = int(2*lag+1)
        lags = np.arange(-lag, lag + 1)

    # get the number of correlation to be computed
    nxcor = _a1xcorPy.a1xcor.get_nxcorr()

    xcorr = np.empty((nxcor, len), dtype='float64')
    ijx = np.empty((nxcor, 2), dtype='int32')


    if (verbose>0):
        print('Computing ',nxcor, 'correlations ')

    if stack is True:
        stack=1
    else:
        stack=0

    if (verbose>=2):
        print('signal length is ',arrayIn.shape[0],'FFT length is ',len)
    xcorr, ijx, ier = _a1xcorPy.a1xcor.compute_xcorr(arrayIn, xcorr, ijx, int(lag), stack, transposed)

    if (verbose>0):
        print('Done')

    if lag_is_second:
        lags = lags * a1.data_header['dt']

    #if first index in array is 0, correct for fortran indexing
    if base ==0:
        ijx -= 1

    return xcorr, lags, ijx, ier

def sort_xcorr(i0, ijx, xcorr=None):
    """
    ##Description
    I, J = sort_xcorr(i0, ijx)
    or
    xcorr_sorted = sort_xcorr(i0, ijx ,xcorr)

    return the xcorr indices I that correspond to correlation couple {i0, ?},
    xcorr[I,:] are the cross-correlation that implies trace i0
    Among these correlation, those with index J need to be time-flipped
    If an xcorr[ncor, nlags] array is given as input, the function returns a sorted version of xcorr

    ##Input:
        src: index of the "source" trace
        ijx: array of indices obtained from the compute_xcorr()

    ##Return
        I, J, [xcorr_sorted]
        I = array of index such that xcorr[I, :] correspond to cross-correlation of couple {src, ?} or {?, src)
        J = array of index xcorr[J,:] where one needs to flip time to get {i0,?} xcorr instead of {?,i0}
        xcorr_sorted[ncor, ntime] = sub-array of xcorr that corresponds only to {i0,?} cross-correlations

    ## Usage example:
        >>#a1 = a <A1Section> instance
        >>xcorr, lags, ijx, ier = xcor.compute_xcorr(a1,lag=lag,stack=True,verbose=0)
        >>
        >># return all the indices {i1} that correspond to correlations {i0,?}
        >>I, J = xcor.sort_xcorr(i0,ijx)
        >>xcorrs = xcor.sort_xcorr(i0,ijx, xcorr)
        >>xcorr_sorted = xcorr[I,:]
        >>xcorr_flipped = numpy.flip(xcorr_sorted[J],axis=1)
        >># xcorrs and xcorr_flipped are identical
    """
    #
    # ijx is an array of trace indices of dimension [ncor x 2]
    # Each row k gives the indices of traces used for the kth correlation
    # That is, for the ncor correlations {i,j}, ijx[:,0] and ijx[:,1] gives the index of traces
    # i and j with i>j

    from numpy import where, unique, flip

    Indx = where(ijx == i0) # list of values {i,j} where i or j is i0
    #Indx[0]
    i1 = Indx[0]    # rows of ijx where the trace index is i0 in col 0 or 1
    i2 = Indx[1]    # cols of ijx where the trace index is i0
    i3 = (i2+1) % 2 # cols of ijx[:,1] where the trace index is not i0, except
                    # if autocorr

    # keep only unique and sorted indices j of correlation {i0,j}
    i1u, i1u_idx = unique(i1, return_index=True)

    # Now find the correlations that needs to be flipped
    # They are those where j is in first column of ijx
    # that is, those where i2 is 1
    i2u = i2[i1u_idx]
    #i3u = i3[i1u_idx]

    # get the list of indices where the xcorr must be flipped in time
    flip_index = where(i2u == 1)
    i1uf = i1u[flip_index]

    if xcorr is not None:
        xcorr_sorted = xcorr[i1u,:]
        xcorr_sorted[i1uf,:]=flip(xcorr_sorted[i1uf,:], axis=1)
        return xcorr_sorted
    else:
        return i1u, i1uf

def get_nxcorr():
    '''
    get the number of xcorr couple of trace that have been registered
    for the computation of cross correlation in the current section
    '''
    from a1das import _a1xcorPy
    return _a1xcorPy.a1xcor.get_nxcorr()

def infos():
    """
    ##Description
    Print information about cross-correlation computations
    """
    from a1das import _a1xcorPy
    _a1xcorPy.a1xcor.display_cxc()




