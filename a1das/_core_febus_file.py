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

# relative tolerance on datation error
# if datation error exceed this tolerance, data are assumed to be unevenly digitized
__time_err_tol__ = 0.01

# DATATION POLICY
# if __true_time_policy__ = True
# Use block start time supplied in Febus file to determine the sample times.
# This implies a non constant time step and sample datation is obtained from the time vector
#
# if __true_time_policy__ = False
# Assume a constant time step given by the Febus header
#
# The policy can be changed dynamically using the a1das.setTrueTimePolicy()
#
__true_time_policy__ = True

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
    from ._a1headers import _A1DataHeader, _A1FileHeader
    from  ._a1das_exception import FileFormatError

    chunk_cache_size = int(1024*1024)
    try:
        f = h5py.File(filename, 'r', rdcc_nbytes = chunk_cache_size)
    except FileNotFoundError:
        print('could not open ', filename, ' as a file')
        raise FileNotFoundError

    chunk_cache_size_is_ok = False
    while not chunk_cache_size_is_ok:
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
        blockTimeNode = f['/' + DASName + '/' + SOURCEName +'/time']
        block_times = blockTimeNode[:]

        # DAS_name = DASName
        prf = float(SrcAttr['PulseRateFreq'] / 1000.)
        try:
            freq_res = float(SrcAttr['FreqRes'] / 1000.)  # A chunk(block) of data is written @ (Freq_res/1000) Hz
        except:
            freq_res = float(SrcAttr['BlockRate'] / 1000.)  # A chunk(block) of data is written @ (BlockRate/1000) Hz
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


        #
        # block structure and length
        # Strain Rate or raw data
        #
        list_name = list(f['/' + DASName + '/' + SOURCEName + '/' + ZONEName])

        dataset_name = list_name[0]

        if dataset_name == 'Strain':
            Name = '/' + DASName + '/' + SOURCEName + '/' + ZONEName + '/' + 'Strain'
            Node = f[Name]
            nb_block = Node.shape[0]
            block_time_size = Node.shape[1]
            block_space_size = Node.shape[2]
            chunk_size = block_time_size * block_space_size
            chunk_byte_size = chunk_size*4
            #CHANNode = None
            type = 'strain'

        elif dataset_name == 'StrainRate':
            Name = '/' + DASName + '/' + SOURCEName + '/' + ZONEName + '/' + 'StrainRate'
            Node = f[Name]
            nb_block = Node.shape[0]
            block_time_size = Node.shape[1]
            block_space_size = Node.shape[2]
            chunk_size = block_time_size * block_space_size
            chunk_byte_size = chunk_size * 4
            #CHANNode = None
            type = 'strain-rate'

        elif dataset_name == 'ch1':
            Name = '/' + DASName + '/' + SOURCEName + '/' + ZONEName + '/' + 'ch'
            Node = [None, None, None, None]
            for i, j in enumerate(range(0, 4)):
                Node[i] = f[Name + str(i + 1)]
            # Assume all phases have same dimensions nblock x space_size x time_size
            nb_block = Node[0].shape[0]
            block_time_size = Node[0].shape[1]
            block_space_size = Node[0].shape[2]
            chunk_size = block_time_size * block_space_size
            chunk_byte_size = chunk_size
            #SRATENode = None
            type = 'raw'

        elif dataset_name == 'Amplitude':
            Name = '/' + DASName + '/' + SOURCEName + '/' + ZONEName + '/' + 'Amplitude'
            Node = f[Name]
            nb_block = Node.shape[0]
            block_time_size = Node.shape[1]
            block_space_size = Node.shape[2]
            chunk_size = block_time_size * block_space_size
            chunk_byte_size = chunk_size * 4
            type = 'amplitude'


        #
        # compare chunk size and chunk_cache_size, if too small, increase chunk_cache_size
        # close file and reopen it with correct cache size
        #
        if chunk_byte_size > chunk_cache_size:
            f.close()
            chunk_cache_size = chunk_byte_size
            f = h5py.File(filename, 'r', rdcc_nbytes=chunk_cache_size)
            print('Febus file: adjusintg cache size to ',chunk_byte_size/1024/1024,' Mbytes')
        else:
            print('Febus file: chunk_size is ',chunk_byte_size,' and cache size ',chunk_cache_size)
            chunk_cache_size_is_ok = True

    #
    # check time continuity
    #
    if len(block_times) > 1:
        # compute time lapse between consecutive data block
        dt_block = np.diff(block_times)
        # average time between consecutive blocks
        dt_mean = np.mean(dt_block) #should be freq_res in theory
        time_errs = np.count_nonzero((dt_block-dt_mean) > __time_err_tol__*dt_mean)
        if time_errs >0:
            print('open_febus_file: WARNING: datation error detected on block time. Block lapse time exceeds tolerance')

    #
    # Block (chunk) info
    #
    block_info = {'block_time_size': block_time_size, 'chunk_time_length_in_sec': block_time_size*dt ,'chunk_space_size': block_space_size, 'chunk_size': chunk_size,
                  'nb_block': nb_block, 'block_times': block_times}

    #
    # distance information
    #
    Extent = ZneAttr['Extent']
    # ZI_start = Origin[0] + Extent[0] * dx
    # ZI_end = Origin[0] + Extent[1] * dx
    nspace = block_space_size
    ospace = Origin[0] + Extent[0] * dx
    #dist=np.linspace(ospace,ospace +nspace*dx, nspace)
    dist = np.arange(0,nspace)*dx

    #
    # time information
    #
    #
    # !!!! At the time we open the file, we don't know yet how we read it.
    # - if we read half block, time starts at block_time_size/4*dt
    # - if we read full block, time starts at 0
    #
    otime = Origin[1]  # WARNING WARNING actually this is not set in the header, timing could be obtained
    otime = block_times[0] #placeholder
    # from it,e.g. Posix time of the recording start

    # number of time samples assuming that we read with skip option = False
    # The true value of npts will be computed later by A1File._get_time_bounds()
    npts = nb_block * int(block_time_size / 2)
    # time assuming constant time step
    # time_cst = np.arange(0,npts)*dt + int(block_time_size / 4)*dt
    # time assuming unevenly time step
    time = _compute_true_time(nb_block, int(block_time_size / 2), dt, block_times) + int(block_time_size / 4)*dt
    time -= otime

    # create _A1FileHeader instance
    a1fh = _A1FileHeader(freq_res, block_info, Node, 'febus', f, filename)

    # create _A1DataHeader instance
    hdr_dict={'gauge_length':gauge_length, 'sampling_res':sampling_res, 'prf':prf, 'derivation_time': derivation_time,
         'axis1':'time', 'axis2':'space', 'dt':dt, 'ntime':npts, 'otime':otime, 'dx':dx, 'nspace':nspace,
         'ospace': ospace, 'data_type': type, 'time':time, 'dist':dist }

    a1dh =  _A1DataHeader (hdr=hdr_dict)

    return a1fh, a1dh


#
# ====================================    READ_FEBUS_FILE()  =======================================
#
def read_febus_file(f, block=None, drange=None, trange=None, ddecim=1, skip=True, verbose=0, float_type='float64'):
    """
    ## Description
    This function is supposed to be called from core.read()

    read_febus_file(f, block=None, trange=None, drange=None,  ddecim=1, skip=True, verbose=0)

    Read all or part of a DAS strain[rate] (not raw) file into memory

    Several options:
    1) read by block number block = [start, end] (end included)
    2) read by time range trange = [start, end]  (end included, in sec)
    3) real all if trange = block = None
    4) read positions comprised between drange = [start, end] or all if drange = None

    ## Input:
    a1 = an A1File class instance<br>
    block = a list of 1 or 2 elements: [one_block_to_read] or [first_block, last_block], or None<br>
    trange = a list of time range to read, see `a1das.core.parse_trange` for format details<br>
    drange = a list of distance range to read: [dmin, dmax], or None<br>
    ddecim = space decimation or 1<br>
    skip = read all blocks or 1 over 2 (faster but read only Even number of blocks)<br>
    float_type= float type on output, data are native float32<br>
    verbose = >= 0<br>

    return:
      data_header (_A1DataHeader class instance), data (numpy 2D ndarray)

    """

    from ._a1das_exception import DataTypeError
    from .core import parse_trange, parse_drange

    hd = f.file_header
    block_times = hd.block_info['block_times']

    if f.data_header['data_type'] == 'raw':
        RAW = True
    else:
        RAW = False

    if float_type == 'float64':
        npftype = np.float64
    elif float_type== 'float32':
        npftype = np.float32
    else:
        raise DataTypeError('Wrong float type '+ float_type )

    dt = f.data_header['dt']

    #
    # parse trange
    #
    trange = parse_trange(f.data_header, trange)
    drange = parse_drange(f.data_header, drange)


    #
    # read according to block or to trange

    # indices for chunk (block) of data
    # indices for first and last time record in the first and last block
    # Vector of times in range[from_time, to_time] with tdecim
    first_block, last_block, step_block, \
        time_out, block_indices = f._get_time_bounds(block=block, trange=trange, skip=skip)

    #
    # ---------------------------  compute distance bounds and indices ----------------
    #

    dist_out, dist_ix, dist_in = f._get_space_bounds(drange, ddecim)

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
    if not RAW:
        data = np.empty((output_time_size, output_space_size), npftype, 'C')
    else:
        data = np.empty((4,output_time_size, output_space_size), np.int8, 'C')
    time = np.empty((output_time_size,))
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

        #read data buffer and sample date
        if RAW:
            for i in range(0,4):
                buff = hd.node[i][block, block_start:block_end, :].astype(float_type)
                if np.isnan(buff).any():
                    np.nan_to_num(buff,copy=False, nan=0.)
                block_time_size = buff.shape[0]
                data[i,offset: offset + block_time_size, :] = buff[:, dist_ix.start:dist_ix.stop:dist_ix.step]
        else:
            buff = hd.node[block, block_start:block_end, :].astype(float_type)
            if np.isnan(buff).any():
                np.nan_to_num(buff,copy=False, nan=0.)
            block_time_size = buff.shape[0]
            # copy in buff_1 decimated space samples and all time samples
            data[offset: offset + block_time_size, :] = buff[:, dist_ix.start:dist_ix.stop:dist_ix.step]


        time_buff = np.arange(block_start, block_end) * dt + block_times[block]
        time[offset: offset + block_time_size] = time_buff
        offset += block_time_size

    # initialize data header section with the file data header values
    data_header = f.data_header.copy()
    # adjust time and space values
    data_header.set_item(nspace=len(dist_out))
    data_header.set_item(dist=dist_out)
    if __true_time_policy__ is True:
        data_header.set_item(ntime=len(time))
        data_header.set_item(otime=time[0])
        data_header.set_item(time=time-time[0])
    else:
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

    buff = hd.node[block, :, :].astype('float32')

    return  buff


def close_das_file(hd):
    """

    """
    hd.fd.close()

def _compute_true_time(nblock, npts, dt, timeb):
    """
    compute the time vector for nblock consecutive block of npts samples.
    Time step between samples is dt, and each block start time is given by timeb[] array
    """
    time = np.ndarray((nblock*npts,))
    brange=np.arange(0,npts)*dt
    for i in range(0,nblock):
        time[i*npts:(i+1)*npts] = timeb[i] + brange

    return time

def _set_true_time_policy(policy=True):
    """
    see definition in a1das.core.set_true_time_policy()
    """
    __true_time_policy__ = policy

def _true_time_policy():
    """
    return the status of the true_time_policy
    """
    return __true_time_policy__