#
#
#     Module core_reducted_file
#
# Everything needed to open, read, close a reducted file, i.e. a file containing das data reducted
# from a Febus file
#
# structure of the file
#
#
# procedures:
#               open_reducted_file()
#               read_reducted_file()
#               close_reducted_file()
#
__version_reducted_format__='1.0'

__compression_pgm__='lzf'

#
# ============================================ open_reducted_file() ============================
#
def open_reducted_file(filename):
    """

    """
    from ._a1headers import _A1DataHeader, _A1FileHeader
    from  ._a1das_exception import FileFormatError
    import h5py

    chunk_cache_size = int(1024 * 1024)
    try:
        f = h5py.File(filename, 'r',rdcc_nbytes = chunk_cache_size)
    except FileNotFoundError:
        print('could not open ',filename,' as a file')
        raise FileNotFoundError

    #
    # read header
    # A reduction file has a group 'header' on top of the tree
    # whose attributes contain the {key,value} pair: {"file_type","reducted_format"}
    #

    try:
        hdr_attr=f['header'].attrs
        file_type = hdr_attr['file_type']
        if file_type != 'reducted_format':
            raise FileFormatError('not a <reducted file> format or old format. Read directly hdf5 or re-run reduction')
    except:
        raise FileFormatError('not a reducted file format, or old format. Read directly hdf5 or re-run reduction')
    version = hdr_attr['version']
    #
    # read all remaining header attributes and push them into a dictionary
    #
    hdr_dict={}
    for key,value in hdr_attr.items():
        if key != 'file_type' and key != 'version':
            hdr_dict[key] = value

    #
    # read axis vectors <time> and <distance>
    #
    hdr_dict['time'] = f['/time'][:]
    hdr_dict['dist'] = f['/distance'][:]

    #
    # manage old format for axis label and transposition
    #
    if 'data_axis' in hdr_dict:
        if hdr_dict['data_axis'] == 'time_x_space':
            hdr_dict['axis1'] = 'time'
            hdr_dict['axis2'] = 'space'
        elif hdr_dict['data_axis'] == 'space_x_time':
            hdr_dict['axis1'] = 'space'
            hdr_dict['axis2'] = 'time'
        del hdr_dict['data_axis']
    if 'transposition' in hdr_dict:
        if hdr_dict['transposition'] == 1:
            hdr_dict['axis1'] = 'space'
            hdr_dict['axis2'] = 'time'
        else:
            hdr_dict['axis1'] = 'time'
            hdr_dict['axis2'] = 'space'
        del hdr_dict['transposition']

    #
    # check if chunk cache size is large enough
    #
    chunk_cache_size_is_ok = False
    while not chunk_cache_size_is_ok:
        # transposed data
        if hdr_dict['axis1'] == 'space':
            try:
                chunk_size = f['Traces/0'].chunks
                chunk_size = chunk_size[0] * 4 #float32 datatype is 4 bytes
            except:
                chunk_size = f['section'].chunks
                chunk_size = chunk_size[0] * chunk_size[1] * 4 #float32 datatype is 4 bytes
        #not transposed data
        else:
            chunk_size = f['section'].chunks
            chunk_size = chunk_size[0] * chunk_size[1] * 4 #float32 datatype is 4 bytes

        if chunk_size > chunk_cache_size:
            f.close()
            chunk_cache_size = chunk_size
            f = h5py.File(filename, 'r', rdcc_nbytes=chunk_cache_size)

        else:
            chunk_cache_size_is_ok = True

    #
    # Fill data and file class headers
    #
    a1dh = _A1DataHeader(hdr = hdr_dict)
    a1fh = _A1FileHeader(None, None, None, 'reducted', f, filename)

    return a1fh, a1dh
#
# ============================================ READ_REDUCTED_FILE() ============================
#
def read_reducted_file(f1, drange=None, trange=None, verbose=None, float_type='float64',tdecim=None,ddecim=None):

    """
    Read a file reducted from original febus file. This function is to be called from a1file.read()

    input:
    - f1 = <A1File> class instance
    - drange = (list/tuple) distance range in meter from 1st sample [start, end]/[value] or (range) index range
                (default = None, read all)
    - trange = (list/tuple) time range, see `core.parse_trange` for details bout format
               (default = None, read all)
    - verbose not use yet
    - float_type = (float type) force type conversion (default=float64)
    - tdecim = (int) time decimation
    - ddecim = (int) space decimation
    """
    from .core import parse_trange, parse_drange
    from numpy import empty, nanargmin, abs
    from scipy import signal
    from numpy import arange
    from  ._a1das_exception import FileFormatError, WrongValueError
    from obspy.core import UTCDateTime

    # H5 file handle
    f = f1.file_header.fd

    if f1.data_header['data_type'] == 'raw':
        raise FileFormatError('reducted raw format not supported for reading. Read directly using h5py')

    time = f1.data_header['time']
    dist = f1.data_header['dist']
    nspace = f1.data_header['nspace']
    ntime = f1.data_header['ntime']
    dt = f1['dt'] # or f1.data_header['dt']
    if tdecim is None:
        tdecim = 1
    if ddecim is None:
        ddecim = 1

    # drange
    drange = parse_drange(f1.data_header, drange)
    if drange is None:
        d1 = 0
        d2 = nspace
    else:
        d1 = nanargmin(abs(dist - drange[0]))
        d2 = nanargmin(abs(dist - drange[1])) #+ 1
        nspace = d2 - d1

    # parse trange and convert to index
    trange = parse_trange(f1.data_header, trange)
    if trange is None:
        t1 = 0
        t2 = ntime
    else:
        t1 = nanargmin(abs(time - trange[0]))
        t2 = nanargmin(abs(time - trange[1])) #+ 1
        ntime = t2 - t1

    #
    # data have been transposed from original
    #
    if f1['axis1'] == 'space':
        #
        # data transposed from original Febus
        #

        #
        # write vebose messages
        #
        if verbose>=1:
            print('reading transposed data')
            try:
                print('chunk size for trace by trace storage',f['Traces/0'].chunks)
            except:
                print('chunk size for array storage',f['section'].chunks)


        if tdecim != 1:
            if tdecim < 0 or tdecim >= ntime / 2:
                raise WrongValueError('wrong value for time decimation')
            f_nyquist = 1. / f1['dt'] / 2.
            f_corner = 0.7 * f_nyquist / tdecim
            sos = signal.butter(6, f_corner / f_nyquist, 'lowpass', output='sos')
            ntime = len(arange(t1,t2,tdecim))
            #tmp = empty((nspace, t2-t1), dtype=float_type)
            section = empty((nspace, ntime), dtype=float_type)
            for i, ix in enumerate(range(d1, d2)):
                try: # 1st format
                    tmp = f['Traces/' + str(ix)][t1:t2].astype(float_type)
                except: #2nd format
                    tmp = f['section'][ix,t1:t2].astype(float_type)
                section[i, :] = signal.sosfiltfilt(sos, tmp)[::tdecim]
        else:
            try:
                section = empty((nspace, ntime), dtype=float_type)
                for i, ix in enumerate(range(d1,d2,ddecim)):
                    section[i, :] = f['Traces/'+str(ix)][t1:t2]
            except:
                section = f['section'][d1:d2:ddecim,t1:t2].astype(float_type)
            tdecim=1

        #
        # data have not been transposed from original
        #
    else:
        if verbose>=1:
            print('reading non transposed data')
            try:
                print('chunk size for array storage',f['strain'].chunks)
            except:
                print('chunk size for array storage',f['section'].chunks)

        if tdecim is not None and tdecim !=1 :
            raise WrongValueError('cannot apply time decimation on section of dimension <time_x_space> ie not transposed')

        try:
            section = f['section'][t1:t2,d1:d2:ddecim].astype(float_type) # format >12-2021
        except:
            section = f['strain'][t1:t2, d1:d2:ddecim].astype(float_type) # format < 12-2021
        tdecim=1

    # initialize data_header
    data_header = f1.data_header.copy()
    # fill time and space constants
    data_header.set_item(nspace= nspace)
    data_header.set_item(ntime = ntime)
    data_header.set_item(dist = dist[d1:d2:ddecim])
    data_header.set_item(time = time[t1:t2:tdecim])
    data_header.set_item(dt = dt*tdecim)


    return data_header, section

# ============================================ SAVE_REDUCTED_FILE() ============================
def _save_reducted_file(a1, filename, chunk_space_size = None):
    """
    save_reducted_file(a1, filename)
    Save the A1Section <a1> into the hdf5 filename <filename> in the reducted format
    """
    import h5py
    from ._a1das_exception import WrongValueError

    try:
        fout=h5py.File(filename,'w')
    except:
        raise WrongValueError('Cannot open filename '+filename)

    dhd = a1.data_header._header
    time = a1['time']
    dist = a1['dist']
    ntime = a1['ntime']
    nspace = a1['nspace']

    #
    # create dataset that contains time vector
    #
    fout.create_dataset('time',data=time,dtype='f8')

    #
    # create dataset that contains distance vector
    #
    fout.create_dataset('distance',data=dist,dtype='f8')

    #
    # create groupe '/header' and fill it
    #
    header_grp = fout.create_group('header')
    # specific to reducted file
    header_grp.attrs['file_type']='reducted_format'
    header_grp.attrs['version']=__version_reducted_format__
    #
    # header values associated to section metadata (data_header)
    for key, value in dhd.items():
        if key != 'time' and key != 'dist':
            header_grp.attrs[key] = value

    # transposed data
    if a1.data_header['axis1'] == 'space':
        chunk_time_size = ntime
        if chunk_space_size is None:
            chunk_space_size = int(nspace/10)
        dset = _create_h5_group_transposed_by_section(fout, nspace, ntime, chunk_time_size, chunk_space_size,
                                                              compression=False)
        dset[:,:] = a1.data[:,:]

    else:
        chunk_time_size = min(8192,int(ntime/10))
        chunk_space_size = nspace
        dset = \
            _create_h5_group_not_transposed(fout, nspace, ntime, chunk_time_size, chunk_space_size, compression=False)
        offset=0
        for i in range(0,ntime,chunk_time_size):
            dset[i:i+chunk_time_size,:]=a1.data[i:i+chunk_time_size,:]
            offset += chunk_time_size
        dset[offset:ntime,:]=a1.data[offset:ntime,:]

    fout.close()

#
# ================================================ _create_h5_group_transposed() =========================
#
def _create_h5_group_transposed_by_trace(fout, fileout, space_size, time_size, chunk_size, compression=True):
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
    chunk_size = h5 time chunck size
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
            dset = grp.create_dataset(str(i), (time_size,),
                                      #chunks=True,
                                      chunks=(chunk_size,),
                                      dtype='f4',compression=__compression_pgm__)
        else:
            dset = grp.create_dataset(str(i), (time_size,),
                                      #chunks=True,
                                      chunks=(chunk_size,),
                                      dtype='f4')

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
    dset = h5py dataset handle to write in
    """

    if compression:
        dset = fout.create_dataset('section', (time_size, space_size),
                                   chunks=(chunk_time_size, chunk_space_size),
                                   #chunks=True,
                                   dtype='f4', compression=__compression_pgm__)
    else:
        dset = fout.create_dataset('section', (time_size, space_size),
                                   chunks=(chunk_time_size, chunk_space_size),
                                   #chunks=True,
                                   dtype='f4')

    return dset


# ================================================ _create_h5_group_transposed() =========================
#
def _create_h5_group_transposed_by_section(fout, space_size, time_size, chunk_time_size, chunk_space_size, compression=True):
    """
    Create the hdf5 group <Section> to write transposed data in a reducted file.
    Data are written in a transposed array

    input:
    fout    = hdf5 file descriptor on output
    space_size = dimension along spatial coordinates
    time_size = dimension along time coordinates
    chunk_time_size, chunk_space_size = h5 chunck size
    compression: True or False

    return:
    dset = h5py dataset handles to write in
    """
    if compression:
        dset = fout.create_dataset('section', (space_size, time_size),
                                   #chunks=True,
                                   chunks=(chunk_space_size, chunk_time_size),
                                   dtype='f4', compression=__compression_pgm__)
    else:
        dset = fout.create_dataset('section', (space_size, time_size),
                                   #chunks=True,
                                   chunks=(chunk_space_size, chunk_time_size),
                                   dtype='f4')


    return dset

