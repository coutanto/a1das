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
#TODO verifier l'etat de l'option ddecim >1 dans les 3 fonctions

__doc__="Functions for converting from Febus format to hdf5 reducted file and perform raw-to-strain[rate] conversion"

#
# ========================================= REDUCTION_TRANSPOSE()  ============================
#
def reduction_transpose(filein, fileout, trange=None, drange=None, tdecim=1, ddecim=1, hpcorner=None, kchunk=10,
                        trace_chunk=None, skip=False, verbose=0, filterOMP=True, no_aliasing=True, compression=False,
                        position=None):
    '''
    ##Description
    Read a DAS section from an hdf5 Febus-Optics format, extract the requested part or all of it,
    perform a transposition of the section, an optional highpass filter in the space domain

    !!! Time decimation is performed with lowpass filtering if requested !!!

    After reduction, data is kept in float32 representation

    ##Input
    filein, fileout: hdf5 (.h5) file to read/write<br>
    trange:  (tuple, list) read time range see `a1das.core.parse_trange` for detail<br>
    default = None, everything is read<br>
    drange:  (tuple, list) distance range in meter (not km) (start,end),see `a1das.core.parse_drange` for detail<br>
    ddecim:  (int) distance decimation factor (default=1)<br>
    tdecim:  (int) time decimation factor (default=1)<br>
    position: (float ndarray[3,:]) 3 x ntraces array with x,y,z position (default=None)<br>
    hpcorner: (float) Corner frequency for High pass spatial filtering (ex. 600m, default=None)<br>
    kchunk:  (int) output_time_chunk_size set to input_time_chunk_size * kchunk (default=10)<br>
    trace_chunk: (int) write trace_chunk traces per space block/chunk (default = ntrace/100)<br>
    skip:    (bool) skip redundant block when reading data (default=False)<br>
    verbose: (int) be verbose (default=0, minimum message level)<br>
    filterOMP: (bool) use the multithread sosfilter binary package instead of scipy.signal.sosfiltfilt (default True)<br>
    no_aliasing: (bool) if tdecim >1, apply a Low pass filter at 0.9 F_nyquist (default True)<br>
    compression: (bool) compress data in H5 files (default False)<br>


    ##Reading transposed reducted output file
    With a1das utilities:

        >>>import a1das
        >>>f=a1das.open(filename,format='reducted') #format is optional
        >>>print(f) # or shortcut f<return>, print file header
        >>>a=f.read() #read all data (or use options) into an a1Section object
        >>>print(a) # print section header parameters
        >>>gl=a['gauge_length'] #get header parameter
        >>>dist=a['dist'] # get distance vector
        >>>time=a['time'] # get time vector
        >>>section=a.data #get 2D float data array
        >>>f.close()

    With h5py:

        >>>import h5py
        >>>f=h5py.File(filename)
        >>>header=f['/header'].attrs   # header dictionnary
        >>>print(header.keys()   )    # list of metadata
        >>>gl=header['gauge_length']  # pick one of metadata
        >>>dist=f['/distance']         # distance vector
        >>>time=f['/time']             # time vector
        >>>section = f['/section'][:,:] # das section
        >>>f.close()


    '''
    import h5py
    from scipy import signal
    from .core import parse_trange, parse_drange
    from ._core_reducted_file import  _create_h5_group_transposed_by_trace, _create_h5_group_transposed_by_section, \
                                    __version_reducted_format__
    from ._core_febus_file import _true_time_policy
    from ._h5util import _get_h5_chunk_cache_size
    from ._a1das_exception import FileFormatError,DataTypeError

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
    f = open(filein, format='febus')
    #
    #  --------------------------- open file for writing
    #
    fout = h5py.File(fileout,'w')


    #
    #  --------------------------- read header ----------------------
    #
    hd = f.file_header
    block_times = hd.block_info['block_times']
    dhd = f.data_header

    if isinstance(hd.node,list):
        raise FileFormatError('Cannot reduce raw data, only strain(rate)')

    #
    # parse (d)trange input format
    #
    drange = parse_drange(f.data_header, drange)
    trange = parse_trange(f.data_header, trange)

    #
    #  --------------------------- block time size
    #  we read the full block and skip redundant blocks
    #  or read half of it and read all blocks
    if skip:
        block_time_size = int(hd.block_info['block_time_size'])
    else:
        block_time_size = int(hd.block_info['block_time_size'] / 2)
    #
    # check whether time decimation factor is ok or not, it must divide
    # the chunk size
    #
    if (block_time_size) % tdecim != 0:
        print('Error: time decimation factor must be a divider of the time chunk size =',block_time_size)
        print('try :')
        for i in range(2,100):
            if (block_time_size) % i == 0:
                print('  tdecim=',i)
        raise DataTypeError('Error: cannot find a correct decimation')


    #
    #   --------------------------- compute time bounds and time vector -----------------
    #

    # indices for chunk (block) of data
    # indices for first and last time record in the first and last block
    # Vector of times in range[from_time, to_time] with tdecim
    first_block, last_block, step_block, \
    time_out, block_indices = f._get_time_bounds(trange=trange, skip=skip, tdecim=tdecim, align_on_block=True)

    #
    #     ---------------------------  compute distance bounds and indices ----------------
    #
    dist_out, dist_ix, dist_in = f._get_space_bounds(drange, ddecim)

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

    if trace_chunk is None:
        out_chunk_space_size = int(f['nspace']/100)
    else:
        out_chunk_space_size = trace_chunk

    print('Original chunk/block is ', block_time_size, 'time samples  x ',input_space_size, ' traces')
    print('Output chunk/block time size is ',out_chunk_time_size,' samples')
    print('       corresponding to ', ncblock, ' original time blocks')
    print('Output space chunk/block size is ',out_chunk_space_size)
    print('output chunk size is ',out_chunk_time_size*out_chunk_space_size*4/1024/1024,' Mbytes')
    #print('chunk cache size: ', _get_h5_chunk_cache_size(f.file_header.fd) / 1024 / 1024,' Mbytes')

    #
    # create groupe '/header' and fill it
    #
    header_grp = fout.create_group('header')
    # header values specific to the h5 reducted file
    header_grp.attrs['file_type']='reducted_format'
    header_grp.attrs['version']=__version_reducted_format__
    # header values associated to section metadata (data_header)
    for key, value in dhd.items():
        if value is not None and key != 'time' and key != 'dist':
            header_grp.attrs[key] = value
    header_grp.attrs['axis1'] = 'space'
    header_grp.attrs['axis2'] = 'time'
    header_grp.attrs['dx'] = dhd['dx']*ddecim
    header_grp.attrs['dt'] = dhd['dt']*tdecim
    header_grp.attrs['nspace'] = output_space_size



    #
    # create dataset that contains distance vector
    #
    fout.create_dataset('distance',data=dist_out,dtype='f8')

    #
    # create groups and datasets
    #
    dset = _create_h5_group_transposed_by_section(fout, output_space_size, output_time_size, out_chunk_time_size,
                                               out_chunk_space_size, compression=compression)

    #
    # --------------------------- loop reading blocks ----------------------
    #
    time = np.empty((output_time_size,))
    buff_out = np.empty((output_space_size, out_chunk_time_size), np.float32, 'C') #modifOC input_space_size=>ouput_space_size
    time_offset = 0 #incremented by output chunk
    time_offset2 = 0 #incremented by input chunk
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

        # copy full current block data into buffer
        # Even if we need only half of the block in the middle
        # reading the full block will allow to do recursive time filtering without side effect
        buff_in = hd.node[block, : , :]
        time_buff = np.arange(block_start, block_end) * dhd['dt'] + block_times[block]

        # highpass space filter if requested
        if hpcorner:
            # replace NaN by 0
            buff_in = np.nan_to_num(buff_in)
            if filterOMP:
                # openMP sosfilter version
                #buff_in[:, :] = _sosfilterPy.sosfiltfilt(sos, buff_in[:, :])
                buff_in = _sosfilterPy.sosfiltfilt_s(sos, buff_in, block_start, block_end, 1)
            else:
                # scipy sequential filter
                #buff_in[:, :] = signal.sosfiltfilt(sos, buff_in[:, :], axis=1)
                buff_in = signal.sosfiltfilt(sos, buff_in, axis=1)
                buff_in[block_start:block_end, :] = signal.sosfiltfilt(sos, buff_in[block_start:block_end, :], axis=1)


        # transpose data, copy is necessary to make sure
        # that the array is transposed in memory
        buff_trans = np.transpose(buff_in[:,dist_ix]).copy() #modifOC buff_in => buff_in[:,dist_ix]

        # increment output buffer block counter
        jchunk = i % ncblock

        if (tdecim>1 and no_aliasing):
            if filterOMP:
                buff_trans = _sosfilterPy.sosfiltfilt_s(sos_time, buff_trans, 0, output_space_size, ddecim)
            else:
                buff_trans[0:output_space_size:ddecim, :] = \
                    signal.sosfiltfilt(sos_time, buff_trans[0:output_space_size:ddecim, :], axis=1)
            #alternative
            #buff_trans[0:output_space_size:ddecim, :] = cusignal.filtering.resample.decimate(buff_trans[0:output_space_size:ddecim, :],tdecim,axis=1)


        # Fill output buffer and perform time decimation
        buff_out[:, jchunk * decim_blk_time_size:(jchunk + 1) * decim_blk_time_size] = \
                 buff_trans[:,block_start:block_start+block_time_size:tdecim]
        time_buff = time_buff[::tdecim]
        l = len(time_buff)
        time[time_offset2:time_offset2+l] = time_buff
        time_offset2 += l

        # when buffer is full , write datasets
        if jchunk == ncblock - 1:
            # output buffer is filled with ncblock time block
            # write it in all spatial dataset
            dset[:, time_offset:time_offset + out_chunk_time_size] = buff_out
            time_offset += out_chunk_time_size

    # end of file reached, write partially filled buffer
    if jchunk != ncblock - 1:
        dset[:, time_offset:time_offset + (jchunk+1)*decim_blk_time_size] = buff_out[:, 0:(jchunk+1)*decim_blk_time_size]


    #
    # create dataset that contains time vector
    #
    if _true_time_policy() is True:
        header_grp.attrs['ntime'] = len(time)
        header_grp.attrs['otime'] = time[0]
        fout.create_dataset('time', data=time-time[0], dtype='f8')
    else:
        header_grp.attrs['ntime'] = output_time_size
        fout.create_dataset('time',data=time_out,dtype='f8')
    f.close()
    fout.close()

#
# ========================================= REDUCTION_NOTRANSPOSE()  ============================
#
def reduction_notranspose(filein, fileout, trange=None, drange=None, tdecim=1, ddecim=1, hpcorner=None, kchunk=10,
                          trace_chunk=None, skip=False, verbose=0, filterOMP=True, compression=False):
    """
    ##Description
    Read a DAS section from an hdf5 Febus-Optics format, extract the requested part or all of it,
    perform an optional highpass filter in the space domain and remove data redundancy

    !!! Time decimation is performed without any lowpass filtering !!!

    After reduction, the records are stored in a single 2D array (time x space) where space = fast axis

    ##Input
    filein, fileout: hdf5 (.h5) file to read/write<br>
    trange:  (tuple, list) read time range see `a1das.core.parse_trange` for detail default = None, everything is read<br>
    drange:  (tuple, list) distance range in meter (not km), see à1das.core.parse_drange` for details,
    (default = None, everything is read)<br>
    ddecim:  (int) distance decimation factor (default=1)<br>
    tdecim:  (int) time decimation factor (default=1)<br>
    hpcorner: (float) Corner frequency for High pass spatial filtering (ex. 600m, default=None)<br>
    kchunk:  (int) the ouput HDF5 chunk size is set to input_time_block_size * kchunk (default=10)<br>
    trace_chunk: (int) output space chunk set to trace_chunk traces (default = ntrace/100)<br>
    skip:    (bool) skip redundant block when reading data (default=False)<br>
    verbose: (int) be verbose (default=0, minimum message level)<br>
    filterOMP: (bool) use the multithread sosfilter binary package instead of scipy.signal.sosfiltfilt (default True)<br>
    compression: (bool) use compression in H5 files (default False)<br>

    ##Reading non transposed reducted output file
    With python using A1File.read(filename,format='reducted')
    or as following example:

        >>>import h5py
        >>>f=h5py.File(filename)
        >>>header=f['header'].attrs   # header dictionnary
        >>>print(header.keys()   )    # list of header keys
        >>>dist=f['/distance']         # distance vector
        >>>time=f['/time']             # time vector
        >>>data=f['/section'][:,:]           # 2D strain[-rate] ndarray [ntime x nspace](float64)

    """

    import h5py
    from .core import open, parse_trange, parse_drange
    from scipy import signal
    from ._core_reducted_file import _create_h5_group_not_transposed,  __version_reducted_format__
    from ._core_febus_file import _true_time_policy
    from ._h5util import _get_h5_chunk_cache_size
    from ._a1das_exception import FileFormatError,DataTypeError

    if filterOMP:
        try:
            from a1das import _sosfilterPy
        except:
            print('could not import sosfilter binary module, switching to scipy')
            filterOMP=False

    #
    #  --------------------------- open file for reading
    #
    f = open(filein, format='febus')
    #
    #  --------------------------- open file for writing
    #
    fout = h5py.File(fileout,'w')


    #
    #  --------------------------- read header ----------------------
    #
    hd = f.file_header
    block_times = hd.block_info['block_times']

    dhd = f.data_header
    if isinstance(hd.node,list):
        raise FileFormatError('Cannot reduce raw data, only strain(rate)')

    #
    # parse (d)trange input format
    #
    drange = parse_drange(f.data_header, drange)
    trange = parse_trange(f.data_header, trange)

    #
    #  --------------------------- block time size
    #  we read the full block and skip redundant blocks
    #  or read half of it and read all blocks
    if skip:
        block_time_size = int(hd.block_info['block_time_size'])
    else:
        block_time_size = int(hd.block_info['block_time_size'] / 2)

    #
    # check whether time decimation factor is ok or not, it must divide
    # the chunk size to avoid complex continuity problem between blocks
    #
    if (block_time_size) % tdecim != 0:
        print('Error: time decimation factor must be a divider of the time chunk size =',block_time_size)
        print('try :')
        for i in range(2,100):
            if (block_time_size) % i == 0:
                print('  tdecim=',i)
        raise DataTypeError('Error: cannot find a correct decimation')

    #
    #   --------------------------- compute time bounds and time vector -----------------
    #

    # indices for chunk (block) of data
    # indices for first and last time record in the first and last block
    # Vector of times in range[from_time, to_time] with tdecim
    first_block, last_block, step_block, \
    time_out, block_indices = f._get_time_bounds(trange=trange, skip=skip, tdecim=tdecim, align_on_block=True)

    #
    #     ---------------------------  compute distance bounds and indices ----------------
    #
    dist_out, dist_ix, dist_in = f._get_space_bounds(drange, ddecim)

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

    out_chunk_time_size = int(block_time_size * kchunk)
    if out_chunk_time_size > output_time_size:
        out_chunk_time_size = output_time_size
        kchunk = int(out_chunk_time_size / block_time_size) + 1
    ncblock = kchunk * tdecim
    tmp_buff_time_size = block_time_size * ncblock

    if trace_chunk is None:
        out_chunk_space_size = int(f['nspace']/100)
    else:
        out_chunk_space_size = trace_chunk

    print('Original chunk/block is ', block_time_size, 'time samples * ',input_space_size, ' traces')
    print('Output chunk/block time size is ',out_chunk_time_size,' samples')
    print('       corresponding to ', ncblock, ' original time blocks')
    print('Output space chunk/block size is ',out_chunk_space_size)
    print('output chunk size is ', out_chunk_time_size * out_chunk_space_size * 4 / 1024/1024, ' Mbytes')
    #print('chunk cache size: ',_get_h5_chunk_cache_size(f.file_header.fd)/1024/1024,' Mbytes')
    #
    # create groupe '/header' and fill it
    #
    header_grp = fout.create_group('header')
    # header values specific to the h5 reducted file
    header_grp.attrs['file_type']='reducted_format'
    header_grp.attrs['version']=__version_reducted_format__
    # header values associated to section metadata (data_header)
    for key, value in dhd.items():
        if value is not None and key != 'time' and key != 'dist':
            header_grp.attrs[key] = value
    header_grp.attrs['axis1']= 'time'
    header_grp.attrs['axis2'] = 'space'
    header_grp.attrs['dx']=dhd['dx']*ddecim
    header_grp.attrs['nspace']=output_space_size
    header_grp.attrs['dt']=dhd['dt']*tdecim
    header_grp.attrs['ntime']=output_time_size

    # create dataset that contains distance vector
    fout.create_dataset('distance',data=dist_out,dtype='f8')


    # create dataset that contains decimated traces
    dset = _create_h5_group_not_transposed(fout, output_space_size, output_time_size,
                                             out_chunk_time_size, out_chunk_space_size, compression)


    #
    # --------------------------- loop reading blocks ----------------------
    #
    # buff_in contain a full chunk of input data (no decimation)
    time = np.empty((output_time_size,))
    chunk_in  = np.empty((block_time_size, input_space_size), np.float32, 'C')
    buff_1    = np.empty((tmp_buff_time_size, output_space_size), np.float32, 'C')
    time_offset = 0
    time_offset2 = 0
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
        chunk_in[:,:] = hd.node[block, block_start:block_end, :]
        time_buff = np.arange(block_start, block_end) * dhd['dt'] + block_times[block]

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
        time_buff = time_buff[::tdecim]
        l = len(time_buff)
        time[time_offset2:time_offset2+l] = time_buff
        time_offset2 += l

        if jchunk == ncblock - 1:
            # write it in all spatial dataset and decimate here for time decimation
            buff_2 = np.squeeze(buff_1[0:tmp_buff_time_size:tdecim, :])
            dset[time_offset:time_offset + out_chunk_time_size, :] = buff_2
            time_offset += out_chunk_time_size

    # end of file reached, write partially filled buffer
    if jchunk != ncblock - 1:
        buff_2 = np.squeeze(buff_1[0:(jchunk + 1) * block_time_size:tdecim, :])
        dset[time_offset:time_offset + (jchunk+1)*decim_blk_time_size, :] = buff_2

    #
    # create dataset that contains time vector
    #
    if _true_time_policy() is True:
        header_grp.attrs['ntime'] = len(time)
        header_grp.attrs['otime'] = time[0]
        fout.create_dataset('time', data=time-time[0], dtype='f8')
    else:
        fout.create_dataset('time',data=time_out,dtype='f8')

    fout.close()
    f.close()


#
# ========================================= RAW2STRAIN()  ============================
#
def raw2strain(filein, fileout, GL, DT, order_time=2, order_space=2, trange=None, drange=None, tdecim=1, ddecim=1,
               kchunk=10, trace_chunk = None, skip=False, verbose=0, compression=False, transpose=True):
    """
    ##Description
    Read a DAS section from an hdf5 Febus-Optics raw format, extract the requested part or all of it,
    remove data redundancy, compute the strain[rate] and store it in a reducted file format

    !!! Decimation are performed without any lowpass filtering !!!

    After reduction, the records are stored in a 2D array [time x space] where space = fast axis

    ##input
    filein, fileout: hdf5 (.h5) file to read/write<br>
    GL: (float) Gauge length (meter)<br>
    DT: (float) Derivation time (second)<br>
    order_time:  (int) finite derivation order in time, no derivation if set to 0 (default 2)<br>
    order_space: (int) finite derivation order in space (default 2)<br>
    trange:  (tuple, list) read time range see `a1das.core.parse_trange` for detail,
     default = None, everything is read<br>
    drange:  (tuple, list) distance range in meter (not km), see à1das.core.parse_drange` for detail,
    (default = None, everything is read)<br>
    ddecim:  (int) read distance decimation factor (default=1)<br>
    tdecim:  (int) read time decimation factor (default=1)<br>
    kchunk:  (int) the ouput HDF5 time chunk size is set to kchunk * input time chunk size (default=10)<br>
    trace_chunk: (int) output space chunk set to trace_chunk traces (default = ntrace/100)<br>
    skip:    (bool) skip redundant block when reading data (default=False)<br>
    verbose: (int) be verbose (default=0, minimum message level)<br>
    compression: (bool) use compression on writing (default False)<br>
    transpose: (bool, default True) write a reducted file that is transposed (space x time) or not transposed (time x space)<br>

    ##Reducted Output format
    Use a minimum HDF5 structure with one group and several Datasets

        group  /
        header attribute = header values
        dataset distance = distance vector
        dataset time =     time vector
        datasets for data depend if data are transposed or not from the original [time x space] ordering

    ##Reading non transposed reducted output file
    With python using A1File.read(filename,format='reducted')
    or as following example:

        >>>import h5py
        >>>f=h5py.File(filename)
        >>>header=f['/header'].attrs   # header dictionnary
        >>>print(header.keys()   )     # list of metadata
        >>>dist=f['/distance']         # distance vector
        >>>time=f['/time']             # time vector
        >>>data=f['/section'][:,:]     # 2D strain[rate] ndarray [ntime x nspace]

    ##Reading transposed reducted output file
    With python using A1File.read(filename,format='reducted')
    or as following example:

        >>>import h5py
        >>>f=h5py.File(filename)
        >>>header=f['/header'].attrs   # header dictionnary
        >>>print(header.keys()   )    # list of metadata
        >>>dist=f['/distance']         # distance vector
        >>>time=f['/time']             # time vector
        >>>section=f['/section'][:,:]  #2D ndarray [nspace x ntime] (float64)

    """

    import h5py
    import numpy as np
    from .core import open, parse_trange, parse_drange
    from a1das import _raw2strPy
    from ._a1das_exception import WrongValueError,DataTypeError
    from ._core_reducted_file import _create_h5_group_not_transposed, _create_h5_group_transposed_by_trace, \
                                     _create_h5_group_transposed_by_section, __version_reducted_format__
    from ._core_febus_file import _true_time_policy

    #
    #  --------------------------- open file for reading
    #
    f = open(filein, format='febus')
    #
    #  --------------------------- open file for writing
    #
    fout = h5py.File(fileout, 'w')

    #
    #  --------------------------- read header ----------------------
    #
    hd = f.file_header
    block_times = hd.block_info['block_times']
    dhd = f.data_header

    #
    # parse (d)trange input format
    #
    drange = parse_drange(f.data_header, drange)
    trange = parse_trange(f.data_header, trange)

    #
    #  --------------------------- block time size
    #  we read the full block and skip redundant blocks
    #  or read half of it and read all blocks
    if skip:
        block_time_size = int(hd.block_info['block_time_size'])
    else:
        block_time_size = int(hd.block_info['block_time_size'] / 2)

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
        raise DataTypeError('Error: cannot find a correct decimation')

    #
    #   --------------------------- compute time bounds and time vector -----------------
    #

    # indices for chunk (block) of data
    # indices for first and last time record in the first and last block
    # Vector of times in range[from_time, to_time] with tdecim
    # (block, trange, skip, align_on_block, tdecim)
    #
    # !! although tdecim may be larger than 1, set it to one since we need the original time vector for
    # doing the raw to strain conversion. We adjust the time vector to the decimation factor later
    first_block, last_block, step_block, \
    time_out, block_indices = f._get_time_bounds(trange=trange, skip=skip, tdecim=1, align_on_block=True)

    #
    #     ---------------------------  compute distance bounds and indices ----------------
    #
    dist_out, dist_ix, dist_in = f._get_space_bounds(drange, ddecim)

    # check for consistancy
    nmin = GL / dhd['dx'] * 10
    if len(dist_out) < nmin:
        print(
            '\n\n!!!!!!!!!! raw2str needs more points along distance to compute strain, increase drange !!!!!!!!!\n\n')
        raise WrongValueError('wrong drange')

    #
    # ---------------------------   print summary  -------------------
    #
    print('> Initialization')
    print('> --------------')
    print('> Data extraction from [', time_out[0], ' - ', time_out[-1], '] sec and from [', dist_out[0], ' - ',
          dist_out[-1], '] m')
    print('> sampling rate are :', dhd['dt'], ' sec and ', dhd['dx'], ' m')

    #
    # --------------------------- size of data to be written ---------------------------
    #
    #
    output_time_size = int(len(time_out)/tdecim)  # number of time samples with decimation tdecim
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

    out_chunk_time_size = (block_time_size * kchunk)
    if out_chunk_time_size > output_time_size:
        out_chunk_time_size = output_time_size
        kchunk = int(out_chunk_time_size / block_time_size) + 1
    ncblock = kchunk * tdecim

    if trace_chunk is None:
        out_chunk_space_size = int(f['nspace']/100)
    else:
        out_chunk_space_size = trace_chunk

    chunk_size = out_chunk_time_size * out_chunk_space_size


    print('Original chunk/block is ', block_time_size, 'time samples * ',input_space_size, ' traces')
    print('Output chunk/block time size is ',out_chunk_time_size,' samples')
    print('       corresponding to ', ncblock, ' original time blocks')
    print('Output chunk/block space size is ',out_chunk_space_size,' traces')
    print('output chunk size is ', chunk_size * 4 / 1024/1024, ' Mbytes')

    #
    # create groupe '/header' and fill it
    #
    header_grp = fout.create_group('header')
    # header values specific to the h5 reducted file
    header_grp.attrs['file_type'] = 'reducted_format'
    header_grp.attrs['version'] = __version_reducted_format__
    # header values associated to section metadata (data_header)
    for key, value in dhd.items():
        if value is not None and key != 'time' and key != 'dist':
            header_grp.attrs[key] = value
    if transpose:
        header_grp.attrs['axis1'] = 'space'
        header_grp.attrs['axis2'] = 'time'
    else:
        header_grp.attrs['axis1'] = 'time'
        header_grp.attrs['axis2'] = 'space'
    if order_time == 0:
        header_grp.attrs['data_type'] = 'strain'
        header_grp.attrs['derivation_time'] = 0.
    else:
        header_grp.attrs['data_type'] = 'strain-rate'
        header_grp.attrs['derivation_time'] = DT
    header_grp.attrs['dx'] = dhd['dx'] * ddecim
    header_grp.attrs['nspace'] = output_space_size
    header_grp.attrs['dt'] = dhd['dt'] * tdecim
    header_grp.attrs['ntime'] = output_time_size

    #
    # test GL and DT parameters
    #
    if GL < dhd['sampling_res'] / 100:
        raise WrongValueError(
            'Gauge length is smaller than sampling resolutionv' + str(GL) + '<' + str(dhd['sampling_res']))
    if DT < (1. / dhd['prf']):
        raise WrongValueError('time derivation is smaller than 1/Pulse_Rate_Frequency')

    # create dataset that contains distance vector
    fout.create_dataset('distance', data=dist_out, dtype='f8')

    # create dataset that contains strain[rate] section
    if transpose:
        strain_dset = _create_h5_group_transposed_by_section(fout, output_space_size, output_time_size,
                                                  out_chunk_time_size, out_chunk_space_size, compression=compression)
    else:
        strain_dset = _create_h5_group_not_transposed(fout, output_space_size, output_time_size,
                                                      out_chunk_time_size, out_chunk_space_size,
                                                      compression=compression)

    #
    # --------------------------- loop reading blocks ----------------------
    #
    # recall ncblock = kchunk*tdecim
    time = np.empty((output_time_size,))
    buff_in = np.empty((block_time_size * ncblock, output_space_size, 4), np.int8, 'C')
    if transpose:
        # strain = np.empty((output_space_size, time_chunk_size), np.float64, 'C')
        strain = np.empty((output_space_size, block_time_size * ncblock), np.float64, 'C')
    else:
        # strain = np.empty((time_chunk_size, output_space_size), np.float64, 'C')
        strain = np.empty((block_time_size * ncblock, output_space_size), np.float64, 'C')

    time_offset = 0
    time_offset2 = 0
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
        phase1_r = hd.node[0][block, block_start:block_end, :]
        phase1_i = hd.node[1][block, block_start:block_end, :]
        # copy decimated spatial data into temporary buffer
        buff_in[jchunk * block_time_size:(jchunk + 1) * block_time_size, :, 0] = \
            phase1_r[:, dist_ix.start:dist_ix.stop:dist_ix.step]
        buff_in[jchunk * block_time_size:(jchunk + 1) * block_time_size, :, 1] = \
            phase1_i[:, dist_ix.start:dist_ix.stop:dist_ix.step]

        # Phase2, same as above
        phase2_r = hd.node[2][block, block_start:block_end, :]
        phase2_i = hd.node[3][block, block_start:block_end, :]
        buff_in[jchunk * block_time_size:(jchunk + 1) * block_time_size, :, 2] = \
            phase2_r[:, dist_ix.start:dist_ix.stop:dist_ix.step]
        buff_in[jchunk * block_time_size:(jchunk + 1) * block_time_size, :, 3] = \
            phase2_i[:, dist_ix.start:dist_ix.stop:dist_ix.step]
        end = timer()
        if verbose >= 2:
            print('elapsed time when reading one block:', end - start)

        time_buff = np.arange(block_start, block_end) * dhd['dt'] + block_times[block]
        time_buff = time_buff[::tdecim]
        l = len(time_buff)
        time[time_offset2:time_offset2+l] = time_buff
        time_offset2 += l
        #
        # We have filled the temporary buffer
        # convert to strain
        #
        if jchunk == ncblock - 1:
            # start = timer()
            ph1_r = np.squeeze(buff_in[0:block_time_size * ncblock, :, 0]).astype('int32')
            ph1_i = np.squeeze(buff_in[0:block_time_size * ncblock, :, 1]).astype('int32')
            ph2_r = np.squeeze(buff_in[0:block_time_size * ncblock, :, 2]).astype('int32')
            ph2_i = np.squeeze(buff_in[0:block_time_size * ncblock, :, 3]).astype('int32')
            # end = timer()
            # if verbose >= 2:
            #    print('elapsed time when decimating time:', end - start)
            strain, GL, DT = _raw2strPy.raw2strpy(strain, ph1_r, ph1_i, ph2_r, ph2_i,
                                                  dist_out, time_out, GL, DT, order_space, order_time, verbose)
            # write to file
            start = timer()
            # array in native ordering [time: space]
            if not transpose:
                strain_dset[time_offset:time_offset + out_chunk_time_size, :] = strain[::tdecim, :]

            # array in transposed ordering [space: time]
            else:
                strain_dset[:, time_offset:time_offset + out_chunk_time_size] = strain[:,::tdecim]

            end = timer()

            if verbose >= 2:
                print('elapsed time when  writing:', end - start)
            time_offset += out_chunk_time_size

    # end of file reached, write partially filled buffer
    if jchunk != ncblock - 1:
        ph1_r = np.squeeze(buff_in[0:(jchunk + 1) * block_time_size, :, 0]).astype('int32')
        ph1_i = np.squeeze(buff_in[0:(jchunk + 1) * block_time_size, :, 1]).astype('int32')
        ph2_r = np.squeeze(buff_in[0:(jchunk + 1) * block_time_size, :, 2]).astype('int32')
        ph2_i = np.squeeze(buff_in[0:(jchunk + 1) * block_time_size, :, 3]).astype('int32')
        strain, GL, DT = _raw2strPy.raw2strpy(strain, ph1_r, ph1_i, ph2_r, ph2_i, dist_out, time_out,
                                              GL, DT, order_space, order_time, verbose)
        if not transpose:
            strain_dset[time_offset:time_offset + (jchunk + 1) * decim_blk_time_size, :] = \
                strain[0:(jchunk + 1) * block_time_size:tdecim, :]

        # array in transposed ordering [space: time]
        else:
            strain_dset[:, time_offset:time_offset + (jchunk + 1) * decim_blk_time_size] = \
                strain[:,0:(jchunk + 1) * block_time_size:tdecim]

    #
    # update header values with value from raw2strpy that have rounded them
    #
    header_grp.attrs['gauge_length'] = GL
    header_grp.attrs['derivation_time'] = DT
    header_grp.attrs['time_derivation_order'] = order_time
    header_grp.attrs['space_derivation_order'] = order_space

    #
    # create dataset that contains time vector
    #
    if _true_time_policy() is True:
        header_grp.attrs['ntime'] = len(time)
        header_grp.attrs['otime'] = time[0]
        fout.create_dataset('time', data=time-time[0], dtype='f8')
    else:
        fout.create_dataset('time',data=time_out,dtype='f8')

    f.close()
    start = timer()
    fout.close()
    end = timer()
    if verbose >= 1:
        print('elapsed time when closing:', end - start)



