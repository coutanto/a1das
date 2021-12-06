#
# Module to perform the reduction of Febus A1 DAS files and the conversion
# from raw to strain
#
# O. Coutant, ISTerre, UGA, 2020 from matlab & python code by Febus
#
#
# Content:
#
#   I) Reduction remove the redundancy of data, perform optional transposition,
#   apply time and/or space decimation, with/without anti_aliasing,
#   selection of time/space ranges.
#   The data are then written to HDF5 files with simple structure
#      redution_transpose()
#      reduction_notranspose()
#      rawreduction_notranspose()
#      rraw2info()
#
#   II) Conversion from raw to strain is applied without transposition,
#   with time and/or space decimation, WITHOUT anti_aliasing,
#   selection of time/space ranges.
#   Strain is written to HDF5 files with simple structure
#     raw2strain()
#     rraw2strain()
#
#     uses raw2strPy.???.so binary module


import numpy as np

__version_reducted_format__='1.0'

__doc__="Functions for converting from Febus format to hdf5 reducted file and perform raw-to-strain[rate] conversion"
#
# ========================================= REDUCTION_TRANSPOSE()  ============================
#
def reduction_transpose(filein, fileout, trange=None, drange=None, tdecim=1, ddecim=1, hpcorner=None, kchunk=10, skip=True, verbose=0, filterOMP=True, no_aliasing=True, compression=True):
    '''
    ##Description
    Read a DAS section from an hdf5 Febus-Optics format, extract the requested part or all of it,
    perform a transposition of the section, an optional highpass filter in the space domain

    !!! Time decimation is performed with lowpass filtering if requested !!!

    After reduction, the individual trace can be accessed directly through direct hdf5 requests

    ##Input
        filein, fileout: hdf5 (.h5) file to read/write
        trange:  (tuple, list) time range in sec (start,end), (default = None, everything is read)
        drange:  (tuple, list) distance range in meter (not km) (start,end), (default = None, everything is read)
        ddecim:  (int) distance decimation factor (default=1)
        tdecim:  (int) time decimation factor (default=1)
        hpcorner: (float) Corner frequency for High pass spatial filtering (ex. 600m, default=None)
        kchunk:  (int) the ouput HDF5 chunk size is set to input_time_block_size * kchunk (default=10)
        skip:    (bool) skip redundant block when reading data (default=true) !!!! READ Doc
        verbose: (int) be verbose (default=0, minimum message level)
        filterOMP: (bool) use the multithread sosfilter binary package instead of scipy.signal.sosfiltfilt (default True)
        no_aliasing: (bool) if tdecim >1, apply a Low pass filter at 0.9 F_nyquist (default True)
        compression: (bool) compress data in H5 files


    ##Reading transposed reducted output file
    With python using A1File.read(filename,format='reducted')
    or as following example:

        >>>import h5py
        >>>f=h5py.File(filename)
        >>>header=f['header'].attrs   # header dictionnary
        >>>print(header.keys()   )    # list of metadata
        >>>dist=f['distance']         # distance vector
        >>>time=f['time']             # time vector
        >>># Two ways to read all or individual traces
        >>># read trace number 9
        >>>trace9 = f['/Traces/9'][:]
        >>># read the full section
        >>># vsection is a virtual HDF5 2D array
        >>>section=f['/vsection'][:,:]  #2D ndarray [nspace x ntime] (float64)
        >>># read only the first 10 traces
        >>>section=f['/vsection'][0:10,:]

    '''
    import h5py

    from scipy import signal
    if filterOMP:
        try:
            from a1das import _sosfilterPy
        except:
            print('could not import sosfilter binary module, switching to scipy')
            filterOMP=False
    from .core import open

    #
    #  --------------------------- open file for reading
    #
    a1 = open(filein, format='febus')
    #
    #  --------------------------- open file for writing
    #
    fout = h5py.File(fileout,'w')


    #
    #  --------------------------- read header ----------------------
    #
    hd = a1.file_header
    dhd = a1.data_header
    if hd.srate_node == None:
        print('Cannot reduce raw data, only strain(rate)')
        exit(-1)
    #
    #  --------------------------- block time size
    #  we read the full block and skip redundant blocks
    #  or read half of it and read all blocks
    if skip:
        block_time_size = int(hd.block_info['time_size'])
    else:
        block_time_size = int(hd.block_info['time_size'] / 2)
    #
    # check whether time decimation factor is ok or not, it must divide
    # the chunk size
    #
    if (block_time_size) % tdecim != 0:
        print('Error: time decimation factor must be a divider of the chunk size =',block_time_size)
        print('try :')
        for i in range(2,100):
            if (block_time_size) % i == 0:
                print('  tdecim=',i)
        exit(0)

    #
    #   --------------------------- compute time bounds and time vector -----------------
    #

    # indices for chunk (block) of data
    # indices for first and last time record in the first and last block
    # Vector of times in range[from_time, to_time] with tdecim
    first_block, last_block, step_block, \
    time_out, block_indices = a1._get_time_bounds(trange=trange, skip=skip, tdecim=tdecim, align_on_block=True)

    #
    #     ---------------------------  compute distance bounds and indices ----------------
    #
    dist_out, dist_ix, dist_in = a1._get_space_bounds(drange, ddecim)

    #
    # ---------------------------   print summary  -------------------
    #
    print(' ')
    print('> Data extraction from [', time_out[0], ' - ', time_out[-1], '] sec and from [', dist_out[0], ' - ',
          dist_out[-1], '] m')
    print('> sampling rate are :',dhd['dt'],' sec and ',dhd['dx'],' m')

    #
    # --------------------------- size of data to be written ---------------------------
    #
    output_time_size = len(time_out)
    output_space_size = len(dist_out)
    input_space_size = len(dist_in)

    #
    # ---------------------------- compute optional filter coefficients
    #
    # In space
    if hpcorner:
        k_nyquist = np.pi/dhd['dx']
        k_corner = 2*np.pi/hpcorner
        sos = signal.butter(3, k_corner / k_nyquist, 'highpass', output='sos')
        #scipy.io.savemat('sos.mat', mdict={'sos': sos})
    # In time
    if tdecim >1:
        f_nyquist  = 1./dhd['dt'] / 2.
        f_corner  = 0.7 * f_nyquist / tdecim
        sos_time = signal.butter(6, f_corner / f_nyquist, 'lowpass', output='sos')

    #
    #   --------------------------- write header dataset structure on output file
    #   Create one dataset per distance
    #
    # Each output block contains decim_blk_time_size time samples
    # we will write decim_blk_time_size * tdecim * kchunk time sample per output chunk
    # A chunk is filled when ncblock are read
    decim_blk_time_size = int(block_time_size/tdecim)

    out_chunk_time_size = block_time_size * kchunk
    if out_chunk_time_size > output_time_size:
        out_chunk_time_size = output_time_size
        kchunk = int(out_chunk_time_size / block_time_size)+1
    ncblock = kchunk * tdecim
    tmp_buff_time_size = block_time_size * ncblock

    print('A block is a HDF5 chunk of data')
    print('Original block has ', block_time_size * input_space_size, ' (time x space) values')
    print('Original block time size is ', block_time_size, ' values')
    print('Decimed block time size is ', decim_blk_time_size, ' values')
    print('total time size is ',output_time_size, ' values')
    print('Output block time-size is ',out_chunk_time_size,' values')
    print('       corresponding to ', ncblock, ' original time blocks')

    #
    # create groupe '/header' and fill it
    #
    header_grp = fout.create_group('header')
    header_grp.attrs['file_type']='reducted_format'
    header_grp.attrs['version']=__version_reducted_format__
    header_grp.attrs['transposition']= 1
    header_grp.attrs['data_type']='strainrate'
    header_grp.attrs['gauge_length']=dhd['gauge_length']
    header_grp.attrs['prf']=dhd['prf']
    header_grp.attrs['sampling_res']=dhd['sampling_res']
    header_grp.attrs['derivation_time']=dhd['derivation_time']
    header_grp.attrs['dx']=dhd['dx']*ddecim
    header_grp.attrs['nspace']=output_space_size
    header_grp.attrs['ospace']=dhd['ospace']
    header_grp.attrs['dt']=dhd['dt']*tdecim
    header_grp.attrs['ntime']=output_time_size
    header_grp.attrs['otime']=dhd['otime']

    #
    # create dataset that contains time vector
    #
    fout.create_dataset('time',data=time_out,dtype='f8')

    #
    # create dataset that contains distance vector
    #
    fout.create_dataset('distance',data=dist_out,dtype='f8')

    #
    # create group '/Traces'
    #
    section_list = _create_h5_group_transposed(fout, fileout, output_space_size, output_time_size, out_chunk_time_size,
                                               compression=compression)

    #
    # --------------------------- loop reading blocks ----------------------
    #
    # buff_in  = np.empty((block_time_size, input_space_size), np.float32, 'C')
    buff_out = np.empty((input_space_size, out_chunk_time_size), np.float32, 'C')
    time_offset = 0
    last_block_read = list(range(first_block, last_block, step_block))[-1]
    for i, block in enumerate(range(first_block, last_block, step_block)):
        if verbose >= 1:
            print('    ' + str(block - first_block + 1) + '/' + str(last_block - first_block) + ' blocks', end ='\r')

        # set block indices to read
        if block == first_block:
            block_start = block_indices[0][0]
            block_end = block_indices[0][1]
        elif block == last_block_read:
            block_start = block_indices[2][0]
            block_end = block_indices[2][1]
        else:
            block_start = block_indices[1][0]
            block_end = block_indices[1][1]

        # copy current block data into buffer
        #buff_in[:, :] = hd.srate_node[block, block_start:block_end, :]
        buff_in = hd.srate_node[block, : , :]

        # highpass space filter if requested
        if hpcorner:
            # replace NaN by 0
            #buff_in[np.where(np.isnan(buff_in))] = 0.
            buff_in[np.where(np.isnan(buff_tmp))] = 0.
            #buff_in = np.nan_to_num(buff_in)
            if filterOMP:
                # openMP sosfilter version
                #buff_in[:, :] = _sosfilterPy.sosfiltfilt(sos, buff_in[:, :])
                buff_in = _sosfilterPy.sosfiltfilt_s(sos, buff_in, block_start, block_end, 1)
            else:
                # scipy sequential filter
                #buff_in[:, :] = signal.sosfiltfilt(sos, buff_in[:, :], axis=1)
                buff_in[block_start:block_end, :] = signal.sosfiltfilt(sos, buff_in[block_start:block_end, :], axis=1)


        # transpose data, copy is necessayr to make sure
        # that the array is transposed in memeory
        #buff_trans = np.transpose(buff_in)
        buff_trans = np.transpose(buff_in).copy()

        # increment output buffer block counter
        jchunk = i % ncblock

        if (tdecim>1 and no_aliasing):
            if filterOMP:
                #print(sos_time.shape,buff_trans.shape,sos_time.dtype,buff_trans.dtype)
                buff_trans = _sosfilterPy.sosfiltfilt_s(sos_time, buff_trans, 0, output_space_size, ddecim)
            else:
                buff_trans[0:output_space_size:ddecim, :] = \
                    signal.sosfiltfilt(sos_time, buff_trans[0:output_space_size:ddecim, :], axis=1)
            #alternative
            #buff_trans[0:output_space_size:ddecim, :] = cusignal.filtering.resample.decimate(buff_trans[0:output_space_size:ddecim, :],tdecim,axis=1)
        # Fill output buffer and perform time decimation
        #buff_out[:, jchunk * decim_blk_time_size:(jchunk+1)*decim_blk_time_size] = buff_trans[:, 0:block_time_size:tdecim]
        buff_out[:, jchunk * decim_blk_time_size:(jchunk + 1) * decim_blk_time_size] = \
                 buff_trans[:,0:block_time_size:tdecim]

        # when buffer is full , write datasets
        if jchunk == ncblock - 1:
            # output buffer is filled with ncblock time block
            # write it in all spatial dataset
            for j, jj in enumerate(range(0,output_space_size,ddecim)):
                dset = section_list[j]
                dset[time_offset:time_offset + out_chunk_time_size] = buff_out[jj, :]
            time_offset += out_chunk_time_size

    # end of file reached, write partially filled buffer
    if jchunk != ncblock - 1:
        for j,jj in enumerate(range(0,output_space_size,ddecim)):
            dset = section_list[j]
            dset[time_offset:time_offset + (jchunk+1)*decim_blk_time_size] = buff_out[jj, 0:(jchunk+1)*decim_blk_time_size]

    a1.close()
    fout.close()

#
# ========================================= REDUCTION_NOTRANSPOSE()  ============================
#
def reduction_notranspose(filein, fileout, trange=None, drange=None, tdecim=1, ddecim=1, hpcorner=None, kchunk=10,
                          skip=True, verbose=0, filterOMP=True, compression=True):
    """
    ##Description
    Read a DAS section from an hdf5 Febus-Optics format, extract the requested part or all of it,
    perform an optional highpass filter in the space domain and remove data redundancy

    !!! Time decimation is performed without any lowpass filtering !!!

    After reduction, the records are stored in a single 2D array (time x space) where space = fast axis

    ##Input
        filein, fileout: hdf5 (.h5) file to read/write
        trange:  (tuple, list) time range in sec (start,end), (default = None, everything is read)
        drange:  (tuple, list) distance range in meter (not km) (start,end), (default = None, everything is read)
        ddecim:  (int) distance decimation factor (default=1)
        tdecim:  (int) time decimation factor (default=1)
        hpcorner: (float) Corner frequency for High pass spatial filtering (ex. 600m, default=None)
        kchunk:  (int) the ouput HDF5 chunk size is set to input_time_block_size * kchunk (default=10)
        skip:    (bool) skip redundant block when reading data (default=true)
        verbose: (int) be verbose (default=0, minimum message level)
        filterOMP: (bool) use the multithread sosfilter binary package instead of scipy.signal.sosfiltfilt (default True)
        compression: (bool) use compression in H5 files (default True)

    ##Reading non transposed reducted output file
    With python using A1File.read(filename,format='reducted')
    or as following example:

        >>>import h5py
        >>>f=h5py.File(filename)
        >>>header=f['header'].attrs   # header dictionnary
        >>>print(header.keys()   )    # list of metadata
        >>>dist=f['distance']         # distance vector
        >>>time=f['time']             # time vector
        >>>data=f['strain']           # 2D strain[-rate] ndarray [ntime x nspace](float64)

    """

    import h5py
    from .core import open
    from scipy import signal
    if filterOMP:
        try:
            from a1das import _sosfilterPy
        except:
            print('could not import sosfilter binary module, switching to scipy')
            filterOMP=False

    #
    #  --------------------------- open file for reading
    #
    a1 = open(filein, format='febus')
    #
    #  --------------------------- open file for writing
    #
    fout = h5py.File(fileout,'w')


    #
    #  --------------------------- read header ----------------------
    #
    hd = a1.file_header
    dhd = a1.data_header
    if hd.srate_node == None:
        print('Cannot reduce raw data, only strain(rate)')
        exit(-1)

    #
    #  --------------------------- block time size
    #  we read the full block and skip redundant blocks
    #  or read half of it and read all blocks
    if skip:
        block_time_size = int(hd.block_info['time_size'])
    else:
        block_time_size = int(hd.block_info['time_size'] / 2)

    #
    # check whether time decimation factor is ok or not, it must divide
    # the chunk size
    #
    if (block_time_size) % tdecim != 0:
        print('Error: time decimation factor must be a divider of the chunk size =',block_time_size)
        print('try :')
        for i in range(2,100):
            if (block_time_size) % i == 0:
                print('  tdecim=',i)
        exit(0)

    #
    #   --------------------------- compute time bounds and time vector -----------------
    #

    # indices for chunk (block) of data
    # indices for first and last time record in the first and last block
    # Vector of times in range[from_time, to_time] with tdecim
    first_block, last_block, step_block, \
    time_out, block_indices = a1._get_time_bounds(trange=trange, skip=skip, tdecim=tdecim, align_on_block=True)

    #
    #     ---------------------------  compute distance bounds and indices ----------------
    #
    dist_out, dist_ix, dist_in = a1._get_space_bounds(drange, ddecim)

    #
    # ---------------------------   print summary  -------------------
    #
    print(' ')
    print('> Data extraction from [', time_out[0], ' - ', time_out[-1], '] sec and from [', dist_out[0], ' - ',
          dist_out[-1], '] m')
    print('> sampling rate are :',dhd['dt'],' sec and ',dhd['dx'],' m')

    #
    # --------------------------- size of data to be written ---------------------------
    #
    output_time_size = len(time_out)  # number of time samples with decimation tdecim
    output_space_size = len(dist_out) # number of space samples with decimation ddecim
    input_space_size = len(dist_in)  # number of space samples without decimation
    #input_time_size = block_time_size
    #
    # ---------------------------- compute optional filter coeffcient
    #
    if hpcorner:
        k_nyquist = np.pi/dhd['dx']
        k_corner = 2*np.pi/hpcorner
        sos = signal.butter(3, k_corner / k_nyquist, 'highpass', output='sos')

    #
    #   --------------------------- write header dataset structure on output file
    #   Create a single dataset for all records, stored with a chunk size
    #   approximately equal to : kchunk * original_chunk_size
    #
    # Each block contains decim_block_time_size time samples
    # we will write decim_block_time_size * tdecim * kchunk time sample per chunk
    # A chunk is filled when ncblocks are read
    decim_blk_time_size = int(block_time_size/tdecim)
    decim_blk_space_size = output_space_size

    out_chunk_space_size = decim_blk_space_size
    out_chunk_time_size = int(block_time_size * kchunk)
    if out_chunk_time_size > output_time_size:
        out_chunk_time_size = output_time_size
        kchunk = int(out_chunk_time_size / block_time_size) + 1
    ncblock = kchunk * tdecim
    tmp_buff_time_size = block_time_size * ncblock

    print('A block is a HDF5 chunk of data')
    print('Original block has ', block_time_size * input_space_size, ' (time x space) values')
    print('Decimed block has ', decim_blk_time_size * output_space_size, ' (time x space) values')
    print('Original time block size is ', block_time_size, ' values')
    print('Decimed time block size is ', decim_blk_time_size, ' values')
    print('total time size is ',output_time_size, ' values')
    print('Output block time size is ', out_chunk_time_size, ' values')
    print('       corresponding to ', ncblock, ' original blocks')
    print('Output block has ',out_chunk_time_size * output_space_size, ' (time x space) values')

    #
    # create groupe '/header' and fill it
    #
    header_grp = fout.create_group('header')
    header_grp.attrs['file_type']='reducted_format'
    header_grp.attrs['version']=__version_reducted_format__
    header_grp.attrs['transposition']= 0
    header_grp.attrs['data_type']='strainrate'
    header_grp.attrs['gauge_length']=dhd['gauge_length']
    header_grp.attrs['prf']=dhd['prf']
    header_grp.attrs['sampling_res']=dhd['sampling_res']
    header_grp.attrs['derivation_time']=dhd['derivation_time']
    header_grp.attrs['dx']=dhd['dx']*ddecim
    header_grp.attrs['nspace']=output_space_size
    header_grp.attrs['ospace']=dhd['ospace']
    header_grp.attrs['dt']=dhd['dt']*tdecim
    header_grp.attrs['ntime']=output_time_size
    header_grp.attrs['otime']=dhd['otime']

    # create dataset that contains time vector
    fout.create_dataset('time',data=time_out,dtype='f8')


    # create dataset that contains distance vector
    fout.create_dataset('distance',data=dist_out,dtype='f8')


    # create dataset that contains decimated traces
    records = _create_h5_group_not_transposed(fout, output_space_size, output_time_size,
                                             out_chunk_time_size, out_chunk_space_size, compression)


    #
    # --------------------------- loop reading blocks ----------------------
    #
    # buff_in contain a full chunk of input data (no decimation)
    chunk_in  = np.empty((block_time_size, input_space_size), np.float32, 'C')
    buff_1    = np.empty((tmp_buff_time_size, output_space_size), np.float32, 'C')
    time_offset=0
    last_block_read = list(range(first_block, last_block, step_block))[-1]
    for i, block in enumerate(range(first_block, last_block, step_block)):
        if verbose >= 1:
            print('    ' + str(block - first_block + 1) + '/' + str(last_block - first_block) + ' blocks', end='\r')

        # set block indices to read
        if block == first_block:
            block_start = block_indices[0][0]
            block_end = block_indices[0][1]
        elif block == last_block_read:
            block_start = block_indices[2][0]
            block_end = block_indices[2][1]
        else:
            block_start = block_indices[1][0]
            block_end = block_indices[1][1]

        # copy current block data into buffer
        chunk_in[:,:] = hd.srate_node[block, block_start:block_end, :]

        # highpass space filter if requested
        if hpcorner:
            # replace NaN by 0
            chunk_in[np.where(np.isnan(chunk_in))] = 0.
            if filterOMP:
                # openMP sosfilter version
                chunk_in[:, :] = _sosfilterPy.sosfiltfilt_s(sos, chunk_in[:, :])
            else:
                # scipy sequential filter
                chunk_in[:, :] = signal.sosfiltfilt(sos, chunk_in[:, :], axis=1)


        # Fill output buffer; when it's full, write datasets
        jchunk = i % ncblock
        # copy in buff_1 decimated space samples and all time samples
        buff_1[jchunk * block_time_size:(jchunk + 1) * block_time_size, :] = chunk_in[:, dist_ix.start:dist_ix.stop:dist_ix.step]
        if jchunk == ncblock - 1:
            # write it in all spatial dataset and decimate here for time decimation
            buff_2 = np.squeeze(buff_1[0:tmp_buff_time_size:tdecim, :])
            records[time_offset:time_offset + out_chunk_time_size, :] = buff_2
            time_offset += out_chunk_time_size

    # end of file reached, write partially filled buffer
    if jchunk != ncblock - 1:
        buff_2 = np.squeeze(buff_1[0:(jchunk + 1) * block_time_size:tdecim, :])
        records[time_offset:time_offset + (jchunk+1)*decim_blk_time_size, :] = buff_2

    a1.close()
    fout.close()

#
# ========================================= RAWREDUCTION_NOTRANSPOSE()  ============================
#
def _rawreduction_notranspose(filein, fileout, trange=None, drange=None, tdecim=1, ddecim=1, kchunk=10,
                          skip=True, verbose=0, use_compression=True):
    """
    filein, fileout: hdf5 (.h5) file to read/write
    trange:  time range in sec (start,end), (default = None, everything is read)
    drange:  distance range in meter (not km) (start,end), (default = None, everything is read)
    ddecim:  distance decimation factor (default=1)
    tdecim:  time decimation factor (default=1)
    kchunk:  the ouput HDF5 chunk size is set to input_time_block_size * kchunk (default=10)
    skip:    skip redundant block when reading data (default=true)
    verbose: be verbose (default=0, minimum message level)

    =======================================================

    Read a DAS raw file from an hdf5 Febus-Optics format, extract the requested part or all of it,
    After reduction, the phasis are stored in 4 single 2D array (time x space) where space = fast axis
    using a one byte binary storage per value

    The reducted file is organized as follow
    /
     group   header
     dataset distance = distance vector
     dataset time =     time vector
     dataset phase0, phase1, phase2, phase3 = real,imag part for 1st and second phase


    =======================================================



    """

    import h5py
    from .core import open
    #
    #  --------------------------- open file for reading
    #
    a1 = open(filein,format='febus')
    #
    #  --------------------------- open file for writing
    #
    fout = h5py.File(fileout, 'w')

    #
    #  --------------------------- read header ----------------------
    #
    hd = a1.file_header
    dhd = a1.data_header
    if hd.chan_node == None:
        print('Cannot read strain(rate)')
        exit(-1)
    #
    #  --------------------------- block time size
    #  we read the full block and skip redundant blocks
    #  or read half of it and read all blocks
    if skip:
        block_time_size = int(hd.block_info['time_size'])
    else:
        block_time_size = int(hd.block_info['time_size'] / 2)


    #
    # check whether time decimation factor is ok or not, it must divide
    # the chunk size
    #
    if (block_time_size) % tdecim != 0:
        print('Error: time decimation factor must be a divider of the chunk size =', block_time_size)
        print('try :')
        for i in range(2, 100):
            if (block_time_size) % i == 0:
                print('  tdecim=', i)
        exit(0)

    #
    #   --------------------------- compute time bounds and time vector -----------------
    #

    # indices for chunk (block) of data
    # indices for first and last time record in the first and last block
    # Vector of times in range[from_time, to_time] with tdecim
    first_block, last_block, step_block, \
    time_out, block_indices = a1._get_time_bounds(trange=trange, skip=skip, tdecim=tdecim)

    #
    #     ---------------------------  compute distance bounds and indices ----------------
    #
    dist_out, dist_ix, dist_in = a1._get_space_bounds__(drange, ddecim)

    #
    # ---------------------------   print summary  -------------------
    #
    print(' ')
    print('> Data extraction from [', time_out[0], ' - ', time_out[-1], '] sec and from [', dist_out[0], ' - ',
          dist_out[-1], '] m')
    print('> sampling rate are :',dhd['dt'],' sec and ',dhd['dx'],' m')

    #
    # --------------------------- size of data to be written ---------------------------
    #
    output_time_size = len(time_out)  # number of time samples with decimation tdecim
    output_space_size = len(dist_out)  # number of space samples with decimation ddecim
    input_space_size = len(dist_in)  # number of space samples without decimation
    # input_time_size = block_time_size


    #
    #   --------------------------- write header dataset structure on output file
    #   Create a single dataset for all records, stored with a chunk size
    #   approximately equal to : kchunk * original_chunk_size
    #
    # Each block contains decim_block_time_size time samples
    # we will write decim_block_time_size * tdecim * kchunk time sample per chunk
    # A chunk is filled when ncblocks are read
    decim_blk_time_size = int(block_time_size / tdecim)
    decim_blk_space_size = output_space_size
    chunk_size0 = (decim_blk_time_size * tdecim * kchunk)
    if chunk_size0 > output_time_size:
        chunk_size0 = output_time_size
        kchunk = int(chunk_size0 / block_time_size) + 1
    chunk_size1 = decim_blk_space_size
    chunk_size = chunk_size0 * chunk_size1
    print('Original block size is ', block_time_size * input_space_size, ' values')
    print('Decimed block size is ', decim_blk_time_size * output_space_size, ' values')
    print('Original time block size is ', block_time_size, ' values')
    print('Decimed time block size is ', decim_blk_time_size, ' values')
    print('total time size is ',output_time_size, ' values')
    print('Chunck size is ', chunk_size, ' values')
    ncblock = kchunk * tdecim
    print('       corresponding to ', ncblock, ' original blocks')

    #
    # create groupe '/header' and fill it
    #
    header_grp = fout.create_group('header')
    header_grp.attrs['file_type']='reducted_format'
    header_grp.attrs['version']=__version_reducted_format__
    header_grp.attrs['transposition']= 0
    header_grp.attrs['data_type']='raw'
    header_grp.attrs['gauge_length']=dhd['gauge_length']
    header_grp.attrs['prf']=dhd['prf']
    header_grp.attrs['sampling_res']=dhd['sampling_res']
    header_grp.attrs['derivation_time']=dhd['derivation_time']
    header_grp.attrs['dx']=dhd['dx']*ddecim
    header_grp.attrs['nspace']=output_space_size
    header_grp.attrs['ospace']=dhd['ospace']
    header_grp.attrs['dt']=dhd['dt']*tdecim
    header_grp.attrs['ntime']=output_time_size
    header_grp.attrs['otime']=dhd['otime']

    # create dataset that contains time vector
    fout.create_dataset('time', data=time_out, dtype='f8')

    # create dataset that contains distance vector
    fout.create_dataset('distance', data=dist_out, dtype='f8')

    # create dataset that contains decimated traces
    # phase1, phase2, phase3,phase4
    section_list = []
    for i in range(0,4):
        if use_compression:
            dset = fout.create_dataset('phase'+str(i), (output_time_size, output_space_size),
                                  chunks=(chunk_size0, chunk_size1),
                                  dtype='i1', compression="lzf")
        else:
            dset = fout.create_dataset('phase'+str(i), (output_time_size, output_space_size),
                                  chunks=(chunk_size0, chunk_size1),
                                  dtype='i1')
        # store it in a list
        section_list.append(dset)
    #
    # --------------------------- loop reading blocks ----------------------
    #
    chan_in = np.empty((block_time_size, input_space_size), np.byte, 'C')
    buff_out = np.empty((chunk_size0*tdecim, output_space_size, 4), np.byte, 'C')
    time_offset = 0
    last_block_read = list(range(first_block, last_block, step_block))[-1]
    for i, block in enumerate(range(first_block, last_block, step_block)):
        if verbose >= 1:
            print('    ' + str(block - first_block + 1) + '/' + str(last_block - first_block) + ' blocks', end='\r')

        # set block indices to read
        if block == first_block:
            block_start = block_indices[0][0]
            block_end = block_indices[0][1]
        elif block == last_block_read:
            block_start = block_indices[2][0]
            block_end = block_indices[2][1]
        else:
            block_start = block_indices[1][0]
            block_end = block_indices[1][1]

        # Fill output buffer; when it's full, write datasets
        jchunk = i % ncblock

        for phase in range(0,4):
        # copy current block data into buffer
            chan_in[:, :] = hd.chan_node[phase][block, block_start:block_end, :]
            buff_out[jchunk * block_time_size:(jchunk + 1) * block_time_size, :, phase] = chan_in[:, dist_ix.start:dist_ix.stop:dist_ix.step]

        if jchunk == ncblock - 1:
            # write it in all spatial dataset and decimate here for spatial decimation
            for phase in range(0,4):
                buf_tmp = np.squeeze(buff_out[0:block_time_size*ncblock:tdecim, :, phase])
                dset = section_list[phase]
                dset[time_offset:time_offset + chunk_size0, :] = buf_tmp  # buff_out[0:block_time_size:tdecim, :]
            time_offset += chunk_size0

    # end of file reached, write partially filled buffer
    if jchunk != ncblock - 1:
        for phase in range(0, 4):
            dset = section_list[phase]
            buf_tmp = np.squeeze(buff_out[0:(jchunk + 1) * block_time_size:tdecim, :, phase])
            dset[time_offset:time_offset + (jchunk + 1) * decim_blk_time_size, :] = buf_tmp
            #buff_out[0:(jchunk + 1) * block_time_size:tdecim, 0:input_space_size:ddecim, phase]

    a1.close()
    fout.close()

#
# ========================================= RAW2STRAIN()  ============================
#
def raw2strain(filein, fileout, GL, DT, order_time=2, order_space=2, trange=None, drange=None, tdecim=1, ddecim=1,
               kchunk=1, skip=True, verbose=0, use_compression=True, transpose=False):
    """
    ##Description
    Read a DAS section from an hdf5 Febus-Optics raw format, extract the requested part or all of it,
    remove data redundancy, compute the strain[rate] and store it in a reducted file format

    !!! Decimation are performed without any lowpass filtering !!!

    After reduction, the records are stored in a 2D array [time x space] where space = fast axis

    ##input
        filein, fileout: hdf5 (.h5) file to read/write
        GL: Gauge length (in meter)
        DT: Derivation time (IN SECOND)
        order_time:  finite derivation order in time, no derivation if set to 0 (default 2)
        order_space: finite derivation order in space (default 2)
        trange:  time range in sec (start,end), (default = None, everything is read)
        drange:  distance range in meter (not km) (start,end), (default = None, everything is read)
        ddecim:  distance decimation factor (default=1)
        tdecim:  time decimation factor (default=1)
        kchunk:  the ouput HDF5 chunk size is set to input_time_block_size * kchunk (default=10)
        skip:    skip redundant block when reading data (default=true)
        verbose: be verbose (default=0, minimum message level)
        use_compression: use lz compression on writing (default True)
        transpose: write a reducted file that is transposed (space x time) or not transposed (time x space)

    ##Reducted Output format
    Use a minimalist format in HDF5 with one group and several Datasets

        group  /
        header attribute = metadata
        dataset distance = distance vector
        dataset time =     time vector
        datasets for data depend if data are transposed or not from the original [time x space] ordering

    ##Reading non transposed reducted output file
    With python using A1File.read(filename,format='reducted')
    or as following example:

        >>>import h5py
        >>>f=h5py.File(filename)
        >>>header=f['header'].attrs   # header dictionnary
        >>>print(header.keys()   )    # list of metadata
        >>>dist=f['distance']         # distance vector
        >>>time=f['time']             # time vector
        >>>data=f['strain']           # 2D strain[rate] ndarray [ntime x nspace]

    ##Reading transposed reducted output file
    With python using A1File.read(filename,format='reducted')
    or as following example:

        >>>import h5py
        >>>f=h5py.File(filename)
        >>>header=f['header'].attrs   # header dictionnary
        >>>print(header.keys()   )    # list of metadata
        >>>dist=f['distance']         # distance vector
        >>>time=f['time']             # time vector
        >>># Two ways to read all or individual traces
        >>># read trace number 9
        >>>trace9 = f['/Traces/9'][:]
        >>># read the full section
        >>># vsection is a virtual HDF5 2D array
        >>>section=f['/vsection'][:,:]  #2D ndarray [nspace x ntime] (float64)
        >>># read only the first 10 traces
        >>>section=f['/vsection'][0:10,:]

    """

    import h5py
    from .core import open
    from a1das import _raw2strPy
    from ._a1das_exception import WrongValueError

    #
    #  --------------------------- open file for reading
    #
    a1 = open(filein, format='febus')
    #
    #  --------------------------- open file for writing
    #
    fout = h5py.File(fileout, 'w')

    #
    #  --------------------------- read header ----------------------
    #
    hd = a1.file_header
    dhd = a1.data_header
    if hd.chan_node == None:
        print('Cannot read raw file')
        exit(-1)
    #
    #  --------------------------- block time size
    #  we read the full block and skip redundant blocks
    #  or read half of it and read all blocks
    if skip:
        block_time_size = int(hd.block_info['time_size'])
    else:
        block_time_size = int(hd.block_info['time_size'] / 2)

    #
    # check whether time decimation factor is ok or not, it must divide
    # the chunk size
    #
    if (block_time_size) % tdecim != 0:
        print('Error: time decimation factor must be a divider of the chunk size =', block_time_size)
        print('try :')
        for i in range(2, 100):
            if (block_time_size) % i == 0:
                print('  tdecim=', i)
        exit(0)

    #
    #   --------------------------- compute time bounds and time vector -----------------
    #

    # indices for chunk (block) of data
    # indices for first and last time record in the first and last block
    # Vector of times in range[from_time, to_time] with tdecim
    # (block, trange, skip, align_on_block, tdecim)
    first_block, last_block, step_block, \
    time_out, block_indices = a1._get_time_bounds(trange=trange, skip=skip, align_on_block=True)


    #
    #     ---------------------------  compute distance bounds and indices ----------------
    #
    dist_out, dist_ix, dist_in = a1._get_space_bounds(drange, ddecim)

    #
    # ---------------------------   print summary  -------------------
    #
    print('> Initialization')
    print('> --------------')
    print('> Data extraction from [', time_out[0], ' - ', time_out[-1], '] sec and from [', dist_out[0], ' - ',
          dist_out[-1], '] m')
    print('> sampling rate are :',dhd['dt'],' sec and ',dhd['dx'],' m')

    #
    # --------------------------- size of data to be written ---------------------------
    #
    output_time_size = len(time_out)  # number of time samples with decimation tdecim
    output_space_size = len(dist_out)  # number of space samples with decimation ddecim
    input_space_size = len(dist_in)  # number of space samples without decimation
    # input_time_size = block_time_size

    #
    #   --------------------------- write header dataset structure on output file
    #   Create a single dataset for all records, stored with a chunk size
    #   approximately equal to : kchunk * original_chunk_size
    #
    # Each block contains decim_block_time_size time samples
    # we will write decim_block_time_size * tdecim * kchunk time sample per chunk
    # A chunk is filled when ncblocks are read
    decim_blk_time_size = int(block_time_size / tdecim)
    decim_blk_space_size = output_space_size
    time_chunk_size = (decim_blk_time_size * tdecim * kchunk)
    if time_chunk_size > output_time_size:
        time_chunk_size = output_time_size
        kchunk = int(time_chunk_size / block_time_size) + 1
    space_chunk_size = decim_blk_space_size
    #TODO verifier que les lignes suivantes sont bonnes, il manque pas ddecim dans chunck_size?
    if not transpose:
        chunk_size = time_chunk_size * space_chunk_size
    else:
        chunk_size = time_chunk_size
    print('Original block size is ', block_time_size * input_space_size, ' values')
    print('Decimed block size is ', decim_blk_time_size * output_space_size, ' values')
    print('Original time block size is ', block_time_size, ' values')
    print('Decimed time block size is ', decim_blk_time_size, ' values')
    print('total time size is ',output_time_size, ' values')
    print('Chunck size is ', chunk_size, ' values')
    ncblock = kchunk * tdecim
    print('       corresponding to ', ncblock, ' original blocks')

    #
    # create groupe '/header' and fill it
    #
    header_grp = fout.create_group('header')
    header_grp.attrs['file_type']='reducted_format'
    header_grp.attrs['version']= __version_reducted_format__
    if transpose:
        header_grp.attrs['transposition']=1
    else:
        header_grp.attrs['transposition']=0
    if order_time == 0:
        header_grp.attrs['data_type']='strain'
        header_grp.attrs['derivation_time'] = 0.
    else:
        header_grp.attrs['data_type']='strain-rate'
        header_grp.attrs['derivation_time'] = DT

    header_grp.attrs['gauge_length']=dhd['gauge_length']
    header_grp.attrs['prf']=dhd['prf']
    header_grp.attrs['sampling_res']=dhd['sampling_res']
    header_grp.attrs['dx']=dhd['dx']*ddecim
    header_grp.attrs['nspace']=output_space_size
    header_grp.attrs['ospace']=dhd['ospace']
    header_grp.attrs['dt']=dhd['dt']*tdecim
    header_grp.attrs['ntime']=output_time_size
    header_grp.attrs['otime']=dhd['otime']
    header_grp.attrs['prf']=dhd['prf']
    header_grp.attrs['sampling_res']=dhd['sampling_res']


    #
    # test GL and DT parameters
    #
    if GL<dhd['sampling_res']/100:
        raise WrongValueError('Gauge length is smaller than sampling resolutionv'+str(GL)+'<'+str(dhd['sampling_res']))
    if DT < (1./dhd['prf']):
        raise WrongValueError('time derivation is smaller than 1/Pulse_Rate_Frequency')

    # create dataset that contains time vector
    fout.create_dataset('time', data=time_out, dtype='f8')

    # create dataset that contains distance vector
    fout.create_dataset('distance', data=dist_out, dtype='f8')

    # create dataset that contains strain[rate] section
    if transpose:
        strain_dset = _create_h5_group_transposed(fout, fileout, output_space_size, output_time_size,
                                                     time_chunk_size, compression=use_compression)
    else:
        strain_dset = _create_h5_group_not_transposed(fout, output_space_size, output_time_size,
                                                 time_chunk_size, space_chunk_size, compression=use_compression)

    #
    # --------------------------- loop reading blocks ----------------------
    #
    # TODO verifier les tdecim suprimés dans les 3 lignes dessous, bizarre car déjà pris en compte au dessus
    buff_out = np.empty((time_chunk_size, output_space_size, 4), np.int8, 'C')
    if transpose:
        strain = np.empty((output_space_size, time_chunk_size), np.float64, 'C')
    else:
        strain = np.empty((time_chunk_size, output_space_size), np.float64, 'C')

    time_offset = 0
    from timeit import default_timer as timer
    last_block_read = list(range(first_block, last_block, step_block))[-1]
    for i, block in enumerate(range(first_block, last_block, step_block)):
        if verbose >= 1:
            print('    ' + str(block - first_block + 1) + '/' + str(last_block - first_block) + ' blocks', end='\r')

        # set block indices to read
        if block == first_block:
            block_start = block_indices[0][0]
            block_end = block_indices[0][1]
        elif block == last_block_read:
            block_start = block_indices[2][0]
            block_end = block_indices[2][1]
        else:
            block_start = block_indices[1][0]
            block_end = block_indices[1][1]

        # Fill output buffer; when it's full, write datasets
        jchunk = i % ncblock

        #
        # Read the full bock and compute angle from real and imaginary part
        # phase 1
        start = timer()
        phase1_r = hd.chan_node[0][block, block_start:block_end, :]
        phase1_i = hd.chan_node[1][block, block_start:block_end, :]
        # copy decimated spatial data into temporary buffer
        buff_out[jchunk * block_time_size:(jchunk + 1) * block_time_size, :, 0] = \
                                           phase1_r[:, dist_ix.start:dist_ix.stop:dist_ix.step]
        buff_out[jchunk * block_time_size:(jchunk + 1) * block_time_size, :, 1] = \
                                           phase1_i[:, dist_ix.start:dist_ix.stop:dist_ix.step]
        
        # Phase2, same as above
        phase2_r = hd.chan_node[2][block, block_start:block_end, :]
        phase2_i = hd.chan_node[3][block, block_start:block_end, :]
        buff_out[jchunk * block_time_size:(jchunk + 1) * block_time_size, :, 2] = \
                                           phase2_r[:, dist_ix.start:dist_ix.stop:dist_ix.step]
        buff_out[jchunk * block_time_size:(jchunk + 1) * block_time_size, :, 3] = \
                                           phase2_i[:, dist_ix.start:dist_ix.stop:dist_ix.step]
        end = timer()
        if verbose >= 1:
            print('elapsed time when reading one block:',end-start)

        #
        # We have filled the temporary buffer
        # extract time decimated part and convert to strain
        #
        if jchunk == ncblock - 1:
            start = timer()
            phase1_r = np.squeeze(buff_out[0:block_time_size * ncblock:tdecim, :, 0]).astype('int32')
            phase1_i = np.squeeze(buff_out[0:block_time_size * ncblock:tdecim, :, 1]).astype('int32')
            phase2_r = np.squeeze(buff_out[0:block_time_size * ncblock:tdecim, :, 2]).astype('int32')
            phase2_i = np.squeeze(buff_out[0:block_time_size * ncblock:tdecim, :, 3]).astype('int32')
            end = timer()
            if verbose >= 1:
                print('elapsed time when decimating time:', end - start)
            strain, GL, DT = _raw2strPy.raw2strpy(strain, phase1_r, phase1_i, phase2_r, phase2_i, dist_out, time_out,
                                GL, DT, order_space, order_time, verbose)
            # write to file
            start = timer()
            if not transpose:
                strain_dset[time_offset:time_offset + time_chunk_size, :] = strain

            else:
                # transposition is performed here in python, faster in fortran?
                for j, jj in enumerate(range(0,output_space_size,ddecim)):
                    dset = strain_dset[j]
                    dset[time_offset:time_offset + time_chunk_size] = strain[jj,:]
            end = timer()

            if verbose >= 1:
                print('elapsed time when  writing:', end - start)
            time_offset += time_chunk_size

    # end of file reached, write partially filled buffer
    if jchunk != ncblock - 1:
        phase1_r = np.squeeze(buff_out[0:(jchunk + 1) * block_time_size:tdecim, :, 0])
        phase1_i = np.squeeze(buff_out[0:(jchunk + 1) * block_time_size:tdecim, :, 1])
        phase2_r = np.squeeze(buff_out[0:(jchunk + 1) * block_time_size:tdecim, :, 2])
        phase2_i = np.squeeze(buff_out[0:(jchunk + 1) * block_time_size:tdecim, :, 3])
        strain, GL, DT = _raw2strPy.raw2strpy(strain, phase1_r, phase1_i, phase2_r, phase2_i, dist_out, time_out,
                            GL, DT, order_space, order_time, verbose)
        if not transpose:
            strain_dset[time_offset:time_offset + (jchunk + 1) * decim_blk_time_size, :] = strain[0:(jchunk + 1) * block_time_size:tdecim, :]

        else:
            for j, jj in enumerate(range(0, output_space_size, ddecim)):
                dset = strain_dset[j]
                dset[time_offset:time_offset + (jchunk + 1) * decim_blk_time_size] = strain[jj, 0:(jchunk + 1) * block_time_size:tdecim]

    #
    # update header values with value from raw2str
    #
    header_grp.attrs['gauge_length'] = GL
    header_grp.attrs['derivation_time'] = DT

    a1.close()
    start = timer()
    fout.close()
    end = timer()
    if verbose >= 1:
        print('elapsed time when closing:', end - start)


#
# ========================================= RRAW2INFO()  ============================
#
def _rraw2info(file):
    """
    Read and print infos from a phase binary file obtained after reduction by rawreduction_notranspose()
    file:

    """
    import h5py
    f = h5py.File(file, 'r')

    # read distance and time vector
    distance_node = f['/distance']
    distance = distance_node[:]
    ndist = distance.shape[0]
    time_node = f['/time']
    time = time_node[:]
    ntime = time.shape[0]

    print('file ',file, 'contains:')
    print('time ',time[0],' to ',time[-1],'; dt =',time[1]-time[0], 'ntime =',ntime)
    print('distance ', distance[0], ' to ', distance[-1], '; dx =', distance[1] - distance[0], 'ndist= ',ndist)

    f.close()


#
# ========================================= RRAW2INFO()  ============================
#
def _rraw2strain(file, gauge_length, delta_t, order_time=2, order_space=2, drange=None, trange=None,verbose=0):
    """
    Read  a phase binary file obtained from Raw data after reduction by <raw_reduction_notranspose> and convert to strain

     file: input file in hdf5 format obtained by <rawreduction_notranspose()>
     gauge_length: distance vector
     delta_t: time vector
     order_time: FD order for time
     order_space: FD order for space
     drange: [min, max] range for distance (m), default None
     trange: [min, max] range for time (sec), default None
    :return:
    """

    import h5py
    from a1das import _raw2strPy

    f = h5py.File(file, 'r')

# read distance and time vector
    distance_node = f['/distance']
    distance = distance_node[:]
    ndist = distance.shape[0]
    time_node = f['/time']
    time = time_node[:]
    ntime = time.shape[0]
    dist_min = distance[0]
    dist_max = distance[-1]
    time_min = time[0]
    time_max = time[-1]
    dx = distance[1] - distance[0]
    dt = time[1] - time[0]

# set ranges
    if not drange:
        drange=(0, ndist)
    else:
        drange[0] = int((drange[0] - dist_min)/dx)
        drange[1] = int((drange[1] - dist_min) / dx)
        if drange[0] > ndist and drange[1] > ndist:
            print('wrong time window, > time_max')
            exit()
        if drange[0] < 0 and drange[1] < 0:
            print('wrong time window < time_min')
            exit()
        if drange[0] < 0:
            drange[0] = 0
        if drange[1] > ndist:
            drange[1] = ndist


    if not trange:
        trange=(0, ntime)
    else:
        trange[0] = int((trange[0] - time_min)/dt)
        trange[1] = int((trange[1] - time_min)/dt)
        if trange[0] > ntime and trange[1] > ntime:
            print('wrong time window, > time_max')
            exit()
        if trange[0] < 0 and trange[1] < 0:
            print('wrong time window < time_min')
            exit()
        if trange[0] < 0:
            trange[0] = 0
        if trange[1] > ntime:
            trange[1] = ntime

    time = time[trange[0]:trange[1]]
    distance = distance[drange[0]:drange[1]]

# read  phase
    p1_node = f['/phase0']
    phase1_r = p1_node[trange[0]:trange[1], drange[0]:drange[1]]
    p2_node = f['/phase1']
    phase1_i = p2_node[ trange[0]:trange[1], drange[0]:drange[1]]
    p3_node = f['/phase2']
    phase2_r = p3_node[trange[0]:trange[1], drange[0]:drange[1]]
    p4_node = f['/phase3']
    phase2_i = p4_node[ trange[0]:trange[1], drange[0]:drange[1]]

#convert to strain
    strain = np.empty(phase1_r.shape,dtype='float64')
    _raw2strPy.raw2strpy(strain, phase1_r, phase1_i, phase2_r, phase2_i, drange, trange,
                        gauge_length, delta_t, order_space, order_time, verbose)


    f.close()

    return distance, time, strain

#
# ================================================ _create_h5_group_transposed() =========================
#
def _create_h5_group_transposed(fout, fileout, space_size, time_size, chunk_size, compression=True):
    """
    Create the hdf5 group <Traces> to write transposed data in a reducted file.
    Data are written in different dataset under <Traces> group = <Traces/0>; <Traces/1>; ...

    Data can also be accessed by reference as the  entire dataset named <section>

    Data can also be accessed through a virtual dataset called <vsection>

    input:
    fout    = hdf5 file descriptor on output
    fileout = filename for output (used to add a virtual dataset)
    space_size = dimension along spatial coordinates
    time_size = dimension along time coordinates
    chunk_size = hdf5 chunck size
    compression: True or False

    return:
    section_list = list of h5py dataset handles to write in (one per trace)
    """
    import h5py

    #
    # create group '/Traces'
    #
    grp = fout.create_group('Traces')

    # create dataset that contains reference (!! version >=2.10.0)
    if hasattr(h5py,'ref_dtype'):
        section_ref = fout.create_dataset('section', (space_size,), dtype=h5py.ref_dtype)
    else:
        ref_dtype = h5py.special_dtype(ref=h5py.Reference)
        section_ref = fout.create_dataset('section', (space_size,), dtype=ref_dtype)

    #
    # create datasets that contains traces in group "Traces"
    #
    section_list = []
    for i in range(0, space_size): #modifOC_ATTENTION, ddecim n'est pas pris en compte dans output_space_size ?
        # define a dataset per spatial location
        if compression:
            dset = grp.create_dataset(str(i), (time_size,), chunks=(chunk_size,), dtype='f4',compression="lzf")
        else:
            dset = grp.create_dataset(str(i), (time_size,), chunks=(chunk_size,), dtype='f4', compression="lzf")

        # store it in a list
        section_list.append(dset)

        # set dataset attribute
        if hasattr(h5py, 'ref_dtype'):
            dset.attrs.create('H5PATH','/Traces/'+str(i))

        # make link between reference and dataset
        section_ref[i] = dset.ref

    #
    # create virtual dataset
    #
    layout = h5py.VirtualLayout(shape=(space_size, time_size), dtype='f4')
    for i in range(0, space_size):
        layout[i]= h5py.VirtualSource(fileout, 'Traces/'+str(i),shape=(time_size,))
    fout.create_virtual_dataset('vsection',layout)

    return section_list

#
# ================================================ _create_h5_group_not_transposed() =========================
#
def _create_h5_group_not_transposed(fout, space_size, time_size, chunk_time_size, chunk_space_size, compression=True):
    """
    Create the hdf5 groups to write not transposed data in a reducted file.

    input:
    fout    = hdf5 file descriptor on output
    space_size = dimension along spatial coordinates
    time_size = dimension along time coordinates
    chunk_time_size = hdf5 chunck size along time dimension
    chunk_space_size = hdf5 chunck size along spatial dimension
    compression = bool, True or False

    return:
    dset = h5py dataset handle to write in (one per trace)
    """

    if compression:
        dset = fout.create_dataset('strain', (time_size, space_size),
                                   chunks=(chunk_time_size, chunk_space_size),
                                   dtype='f4', compression="lzf")
    else:
        dset = fout.create_dataset('strain', (time_size, space_size),
                                   chunks=(chunk_time_size, chunk_space_size),
                                   dtype='f4')

    return dset

