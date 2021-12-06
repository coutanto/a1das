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

def open_reducted_file(filename):
    """

    """
    from .core import _A1DataHeader, _A1FileHeader
    from  ._a1das_exception import FileFormatError
    import h5py
    from numpy import arange


    try:
        f = h5py.File(filename, 'r')
    except FileNotFoundError:
        print('could not open ',filename,' as a file')
        raise FileNotFoundError

    #
    # read header
    # A reduction file has a group 'header' on top of the tree

    try:
        hda=f['header'].attrs
        file_type = hda['file_type']
        if file_type != 'reducted_format':
            raise FileFormatError('not a reducted file format')
    except:
        raise FileFormatError('not a reducted file format, or old format. Read directly hdf5 or re-run reduction')

    transpose = hda['transposition']
    if transpose == 0:
        axis = 'time_x_space'
    else:
        axis = 'space_x_time'
    data_type = hda['data_type']
    gauge_length = hda['gauge_length']
    prf = hda['prf']
    sampling_res = hda['sampling_res']
    derivation_time = hda['derivation_time']

    dx = hda['dx']
    nspace=hda['nspace']
    ospace=hda['ospace']
    dt= hda['dt']
    ntime = hda['ntime']
    otime = hda['otime']

    time = f['/time'][:]
    dist = f['/distance'][:]
    #time = arange(0,ntime)*dt+otime
    #dist = arange(0,nspace)*dx+ospace


    #
    # Fill data and file class headers in the DAS class general header
    #
    a1dh = _A1DataHeader(gauge_length=gauge_length, sampling_res=sampling_res, prf=prf, derivation_time=derivation_time,
                           data_type=data_type, data_axis=axis, dt=dt, ntime=ntime, otime=otime, dx=dx, nspace=nspace,
                           ospace=ospace, time=time, dist = dist)
    a1fh = _A1FileHeader(None, None, None, None, 'reducted', f, filename)

    return a1fh, a1dh
#
# ============================================ READ_REDUCTED_FILE() ============================
def read_reducted_file(f1, drange=None, trange=None, verbose=None, float_type='float64'):

    """
    Read a file reducted from original febus file

    input:
    - f1 = <A1File> class instance
    - drange = distance range drange=[start, end] in meters, or None
    - trange = time range trange=[start, end] in sec, or None
    - verbose not use yet
    - float_type = force type conversion (default=float64)
    """
    from numpy import empty, nanargmin, abs
    from  ._a1das_exception import FileFormatError, WrongValueError

    f = f1.file_header.fd

    if f1.data_header['data_type'] == 'raw':
        raise FileFormatError('reducted raw format not supported for reading. Read directly using h5py or rraw2strain')



    time = f1.data_header['time']
    dist = f1.data_header['dist']
    nspace = f1.data_header['nspace']
    ntime = f1.data_header['ntime']

    # drange
    if drange is None:
        d1 = 0
        d2 = nspace
    else:
        d1 = nanargmin(abs(dist - drange[0]))
        if len(drange) == 1:
            d2 = d1 + 1
        else:
            d2 = nanargmin(abs(dist - drange[1])) + 1
        nspace = d2 - d1

    # trange
    if trange is None:
        t1=0
        t2=ntime
    else:
        t1 = nanargmin(abs(time - trange[0]))
        if len(trange) == 1:
            t2 = t1 + 1
        else:
            t2 = nanargmin(abs(time - trange[1])) + 1
        ntime = t2 - t1
    #
    # data have been transposed from original
    #
    if f1.data_header['data_axis'] == 'space_x_time':
        #
        # data transposed from original
        #
        section = empty((nspace, ntime), dtype=float_type)
        for i, ix in enumerate(range(d1,d2)):
            section[i, :] = f['Traces/'+str(ix)][t1:t2]

        #
        # data have not been transposed from original
        #
    else:
        if f1.data_header['data_axis'] != 'time_x_space':
            raise WrongValueError('data_axis in data head er is neither <space_x_time> nor <time_x_space>')
        section = f['strain'][t1:t2,d1:d2].astype(float_type)

    # initialize data_header
    data_header = f1.data_header.copy()
    # fill time and space constants
    data_header.set_item(nspace= nspace)
    data_header.set_item(ntime = ntime)
    data_header.set_item(dist = dist[d1:d2])
    data_header.set_item(time = time[t1:t2])

    return data_header, section

# ============================================ SAVE_REDUCTED_FILE() ============================
def _save_reducted_file(a1, filename):
    """
    save_reducted_file(a1, filename)
    Save the A1Section <a1> into the hdf5 filename <filename> in the reducted format
    """
    import h5py
    from ._a1das_exception import WrongValueError
    from .reduction import _create_h5_group_not_transposed, _create_h5_group_transposed, __version_reducted_format__
    from .core import _A1DataHeader

    try:
        fout=h5py.File(filename,'w')
    except:
        raise WrongValueError('Cannot open filename '+filename)

    nspace = a1['nspace']
    ntime = a1['ntime']
    time = a1['time']
    dist = a1['dist']


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
    header_grp.attrs['file_type']='reducted_format'
    header_grp.attrs['version']=__version_reducted_format__
    header_grp.attrs['data_type']=a1['data_type']
    header_grp.attrs['gauge_length']=a1['gauge_length']
    header_grp.attrs['prf']=a1['prf']
    header_grp.attrs['sampling_res']=a1['sampling_res']
    header_grp.attrs['derivation_time']=a1['derivation_time']
    header_grp.attrs['dx']=a1['dx']
    header_grp.attrs['nspace']=nspace
    header_grp.attrs['ospace']=a1['ospace']
    header_grp.attrs['dt']=a1['dt']
    header_grp.attrs['ntime']=ntime
    header_grp.attrs['otime']=a1['otime']

    # transposed data
    if a1.data_header['data_axis'] == 'space_x_time':
        header_grp.attrs['transposition'] = 1
        chunk_size = ntime
        section_list = _create_h5_group_transposed(fout, filename, nspace, ntime, chunk_size, compression=True)
        for i in range(0,nspace):
            dset = section_list[i]
            dset[:] = a1.data[i,:]

    else:
        header_grp.attrs['transposition'] = 0
        chunk_time_size = min(8192,int(ntime/10))
        chunk_space_size = nspace
        dset = \
            _create_h5_group_not_transposed(fout, nspace, ntime, chunk_time_size, chunk_space_size, compression=True)
        offset=0
        for i in range(0,ntime,chunk_time_size):
            dset[i:i+chunk_time_size,:]=a1.data[i:i+chunk_time_size,:]
            offset += chunk_time_size
        dset[offset:ntime,:]=a1.data[offset:ntime,:]

    fout.close()
