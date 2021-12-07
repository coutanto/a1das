
def sort_xcorr(i0, ijx, xcorr=None):
    """
    I, J [,xcorr_sorted] = sort_xcorr(i0, ijx [,xcorr])
    Get the xcorr indices I that correspond to correlation couple {i0, ?}
    That is xcorr(I,:) are the cross-correlation that implies trace i0
    The array of indices J
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
        >># return all the indices {i1} that correspond to correlations {i0,i1}
        >>I, J, xcorrs =xcor.target_from_src(i0,ijx, xcorr)
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
        return i1u, i1uf, xcorr_sorted
    else:
        return i1u, i1uf