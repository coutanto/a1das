#
# Module core_febus_file
#
# Everything needed to open, read, close a Febus hdf5 strain file and return data to python
# Allows for full reading or reading half data (skip redundancy), space decimation, select time range and space range
#
# O. Coutant, ISTerre, UGA, 2020 from matlab & python code by Febus
#
#
# Content:
#
# Functions:
#               open_febus_file()
#               read_febus_file()
#               read_febus_file_block()
#               close_febus_file()
#
import numpy as np


#
# ====================================    OPEN_FEBUS_FILE()  =======================================
#
def open_febus_file(filename):
    """
    This function is supposed to be called from core.open()
    Open a febus file, fill in some _A1FileHeader and _A1DataHeader
    class informations and return class instances
    input:
        filename = a febus file name
    return:
        _A1FileHeader, _A1DataHeader
     """

    import h5py
    from .core import _A1DataHeader, _A1FileHeader
    from  ._a1das_exception import FileFormatError

    try:
        f = h5py.File(filename, 'r')
    except FileNotFoundError:
        print('could not open ', filename, ' as a file')
        raise FileNotFoundError

    #
    # Organisation:    / DAS / SOURCE / ZONE
    #
    Zone = 1

    # Initialisation
    # Recuperation des donnees d'enetete, dans les attributs du groupe Source1
    #
    DAS = list(f.keys())  # only 1 elemts in this dict.
    try:
        DASName = DAS[0]  # 1st key gives name of groupe
        DASNode = f['/' + DASName]
    except:
        raise FileFormatError('File not in <Febus_H5_format> or wrong versions')

    SOURCES = list(DASNode.keys())
    SOURCEName = SOURCES[0]
    SOURCENode = f['/' + DASName + '/' + SOURCEName]
    SrcAttr = SOURCENode.attrs  # it's a dict.
    ZONEName = 'Zone' + str(Zone)
    ZONENode = f['/' + DASName + '/' + SOURCEName + '/' + ZONEName]
    ZneAttr = ZONENode.attrs

    # DAS_name = DASName
    prf = SrcAttr['PulseRateFreq'] / 1000.
    freq_res = float(SrcAttr['FreqRes'] / 1000.)  # A chunk(block) of data is written @ (Freq_res/1000) Hz
    # Acquisition_length = SrcAttr['FiberLength']
    sampling_res = float(SrcAttr['SamplingRes'])

    try:
        derivation_time = float(ZneAttr['DerivationTime'])
    except KeyError:
        derivation_time = None

    #
    # Space and time origin
    #
    Origin = ZneAttr['Origin']

    #
    # space and time steps
    #
    Spacing = ZneAttr['Spacing']
    dt = Spacing[1] * 1.e-3  # time step in sec (given originally in msec)
    dx = Spacing[0]  # distance step in m

    try:
        gauge_length = float(ZneAttr['GaugeLength'])
    except KeyError:
        gauge_length = 0.

    TIMEName = 'time'
    # TIMENode = f['/' + DASName + '/' + SOURCEName + '/' + TIMEName]

    #
    # block structure and length
    # Strain Rate or raw data
    #
    list_name = list(f['/' + DASName + '/' + SOURCEName + '/' + ZONEName])

    dataset_name = list_name[0]

    if dataset_name == 'Strain':
        SRATEName = '/' + DASName + '/' + SOURCEName + '/' + ZONEName + '/' + 'Strain'
        SRATENode = f[SRATEName]
        nb_block = SRATENode.shape[0]
        block_time_size = SRATENode.shape[1]
        block_space_size = SRATENode.shape[2]
        chunk_size = block_time_size * block_space_size
        CHANNode = None
        type = 'strain'

    elif dataset_name == 'StrainRate':
        SRATEName = '/' + DASName + '/' + SOURCEName + '/' + ZONEName + '/' + 'StrainRate'
        SRATENode = f[SRATEName]
        nb_block = SRATENode.shape[0]
        block_time_size = SRATENode.shape[1]
        block_space_size = SRATENode.shape[2]
        chunk_size = block_time_size * block_space_size
        CHANNode = None
        type = 'strain-rate'

    elif dataset_name == 'ch1':
        CHANName = '/' + DASName + '/' + SOURCEName + '/' + ZONEName + '/' + 'ch'
        CHANNode = [None, None, None, None]
        for i, j in enumerate(range(0, 4)):
            CHANNode[i] = f[CHANName + str(i + 1)]
        # Assume all phases have same dimensions nblock x space_size x time_size
        nb_block = CHANNode[0].shape[0]
        block_time_size = CHANNode[0].shape[1]
        block_space_size = CHANNode[0].shape[2]
        chunk_size = block_time_size * block_space_size
        SRATENode = None
        type = 'raw'

    #
    # Block (chunk) info
    #
    block_info = {'time_size': block_time_size, 'space_size': block_space_size, 'chunk_size': chunk_size,
                  'nb_block': nb_block}

    #
    # distance information
    #
    Extent = ZneAttr['Extent']
    # ZI_start = Origin[0] + Extent[0] * dx
    # ZI_end = Origin[0] + Extent[1] * dx
    nspace = block_space_size
    ospace = Origin[0] + Extent[0] * dx
    #dist=np.linspace(ospace,ospace +nspace*dx, nspace)
    dist = np.arange(0,nspace)*dx + ospace

    #
    # time information
    #
    otime = Origin[1]  # WARNING WARNING actually this is not set in the header, timing could be obtained
    # from it,e.g. Posix time of the recording start
    # number of time samples assuming that we read with skip option = False
    # The true value of npts will be computed later by A1File._get_time_bounds()
    npts = nb_block * int(block_time_size / 2)
    #time_info = {'step': dt, 'npts': npts, 'start': otime,
    #             'end': npts * dt + otime}
    #time = np.linspace(otime, npts * dt + otime, npts)
    time = np.arange(0,npts)*dt + otime

    # create _A1FileHeader instance
    a1fh = _A1FileHeader(freq_res, block_info, SRATENode, CHANNode, 'febus', f, filename)

    # create _A1DataHeader instance
    a1dh = _A1DataHeader(gauge_length=gauge_length, sampling_res=sampling_res, prf=prf, derivation_time=derivation_time,
                           type=type, data_axis='time_x_space',dt=dt, time=time, ntime=npts, otime=otime, dx=dx, dist=dist,
                           nspace=nspace, ospace=ospace, data_type=type)

    return a1fh, a1dh


#
# ====================================    READ_FEBUS_FILE()  =======================================
#
def read_febus_file(a1, block=None, trange=None, drange=None,  ddecim=1, skip=True, verbose=0, float_type='float64'):
    """
    This function is supposed to be called from core.read()

    read_febus_file(a1, block=None, trange=None, drange=None,  ddecim=1, skip=True, verbose=0)

    Read all or part of a DAS strain[rate] (not raw) file into memory

    Several options:
    1) read by block number block = [start, end] (end included)
    2) read by time range trange = [start, end]  (end included, in sec)
    3) real all if trange = block = None
    4) read positions comprised between drange = [start, end] or all if drange = None

    input:
      a1 = an A1File class instance
      block = a list of 1 or 2 elements: [one_block_to_read] or [first_block, last_block], or None
      trange = a list of time range to read: [tmin, tmax], or None
      drange = a list of distance range to read: [dmin, dmax], or None
      ddecim = space decimation or 1
      skip = read all blocks or 1 over 2 (faster but read only Even number of blocks)
      float_type= float type on output, data are native float32
      verbose = >= 0

    return:
      _A1DataHeader (class instance), data (numpy 2D ndarray)

    """

    from ._a1das_exception import DataTypeError

    if a1.data_header['data_type'] == 'raw':
        raise DataTypeError('can only read strain[rate] data')

    if float_type == 'float64':
        npftype = np.float64
    elif float_type== 'float32':
        npftype = np.float32
    else:
        raise DataTypeError('Wrong float type '+ float_type )

    hd = a1.file_header
    #
    # read according to block or to trange

    # indices for chunk (block) of data
    # indices for first and last time record in the first and last block
    # Vector of times in range[from_time, to_time] with tdecim
    first_block, last_block, step_block, \
        time_out, block_indices = a1._get_time_bounds(block=block, trange=trange, skip=skip)

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
    #
    # --------------------------- initialize data ndarray ---------------------------
    #
    output_time_size = len(time_out)
    output_space_size = len(dist_out)

    #
    # --------------------------- loop reading blocks ----------------------
    #
    data = np.empty((output_time_size, output_space_size), npftype, 'C')
    offset = 0
    last_block_read = list(range(first_block, last_block, step_block))[-1]
    for i, block in enumerate(range(first_block, last_block, step_block)):
        if verbose >= 1:
            print('    ' + str(block - first_block + 1) + '/' + str(last_block - first_block) + ' blocks', end='\r')

        # copy current block data into buffer
        if block == first_block:
            block_start = block_indices[0][0]
            block_end = block_indices[0][1]
        elif block == last_block_read:
            block_start = block_indices[2][0]
            block_end = block_indices[2][1]
        else:
            block_start = block_indices[1][0]
            block_end = block_indices[1][1]

        buff = hd.srate_node[block, block_start:block_end, :].astype(float_type)
        if np.isnan(buff).any():
            np.nan_to_num(buff,copy=False, nan=0.)
        block_time_size = buff.shape[0]

        # copy in buff_1 decimated space samples and all time samples
        data[offset: offset + block_time_size, :] = buff[:, dist_ix.start:dist_ix.stop:dist_ix.step]
        offset += block_time_size

    # initialize data header section with the file data header values
    data_header = a1.data_header.copy()
    # adjust time and space values
    data_header.set_item(nspace=len(dist_out))
    data_header.set_item(dist=dist_out)
    data_header.set_item(ntime=len(time_out))
    data_header.set_item(time=time_out)

    return data_header, data

#
# ====================================    __READ_FEBUS_FILE_block()  =======================================
#
def _read_febus_file_block(a1, block=None):
    """
    __read_febus_file_block(A1File, block=None)

    Read a single block (or chunck) of data from a febus file opened by core.open()
    and return it AS FLOAT 32
    THIS IS A PRIVATE FUNCTION

    input:
        a1 = an A1File instance
        block = block index to read
    return:
        a numpy 2D ndarray for float32 containing the data block
    """

    from ._a1das_exception import DataTypeError, WrongValueError

    if block is None:
        raise WrongValueError('read_febus_file_block() requires a block number')

    if a1.data_header['data_type'] == 'raw':
        raise DataTypeError('can only read strain[rate] data')

    hd = a1.file_header

    buff = hd.srate_node[block, :, :].astype('float32')

    return  buff


def close_das_file(hd):
    """

    """
    hd.fd.close()
