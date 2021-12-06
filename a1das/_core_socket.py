#
#
#     Module a1das_socket
#
# Everything needed to open, read, close a Febus socket stream
#
# O. Coutant, ISTerre, UGA, 2021 based on python code by G. Calbris, Febus
#
#
#
#

def open_das_socket(socket_name):
    """
    Open a ZMQ socket to read real-time data flux from a Febus DAS
    input:
        socket_name = a name of type "tcp://adress:port"
    return:
        _A1SocketHeader, _A1DataHeader
    """
    import zmq
    from .core import _A1SocketHeader, _A1DataHeader
    from ._a1das_exception import ZMQSocketError, DataTypeError
    import numpy as np
    import array
    from sys import platform

    # windows/linux settings
    if platform == 'linux' or platform == 'darwin':
        idata_type = 'i'
    else:
        idata_type = 'l'

    # ZMQ Connection
    try:
        context = zmq.Context()
        socket = context.socket(zmq.REQ)
        socket.connect(socket_name)
    except:
        raise ZMQSocketError("Cannot open and/or connect to "+socket_name+" socket")

    #
    # read a first block to get header information
    #

    # send ACQ
    try:
        socket.send(np.double(0))
    except:
        raise ZMQSocketError("cannot send ACQ to "+socket_name+" socket")

    # receive header and data for 1st block
    try:
        message1 = socket.recv()
        #print('message1',len(message1))
        message2 = socket.recv()
        #print('message2', len(message2))
        message3 = socket.recv()
        #print('message3', len(message3))
    except:
        raise ZMQSocketError("cannot receive data from " + socket_name + " socket")

    # Get information
    DATA1 = array.array(idata_type, message1[0:4])
    count = DATA1[0]
    DATA1 = array.array('d', message1[4::])
    otime = DATA1[0]
    time_stamp = DATA1[0]
    #print('open time_stamp',time_stamp)
    DATA2 = array.array('d',message2[0:48])
    dx = DATA2[0]
    dt = DATA2[1] / 1000.
    ospace = DATA2[3]
    DATA2 = array.array(idata_type, message2[48:64])
    nspace = DATA2[1]-DATA2[0] + 1
    ntime = DATA2[3]-DATA2[2] + 1

    dist = (np.arange(0, nspace) * dx + ospace)
    time = np.arange(int(ntime/4),3*int(ntime/4)) * dt
    used_ntime = int(ntime/2)

    if count==4:
        type = 'raw'
    else:
        if ntime>0:
            type = 'strainrate'
        else:
            type = 'energyband'

    if type != 'strainrate':
        raise DataTypeError('socket only support strain rate type')

    dhd = _A1DataHeader(data_axis="time_x_space", type=type, dt=dt, ntime=used_ntime, otime=otime, dx=dx, nspace = nspace,
                                     ospace=ospace, dist=dist, time=time)
    shd = _A1SocketHeader(socket=socket,  socket_name=socket_name, time_stamp=time_stamp, running_time=0.)

    return shd, dhd

def read_das_socket(a1, block=None, drange=None, idrange=None, ddecim=1):
    """

    Read the data from the DAS socket following ZMQ protocol
    Data are read by block(s), the number of block read and stored
    in the data buffer is specified by the <block> argument
    The spatial range of data can be changed using the drange
    and ddecim parameter.

    input:
        block = None, [], [1] => read a single block
                [k,l] => read l-k block
                [k]   => read k blocks
                default is None, read 1 block
        drange = in meter, refer to distance along fiber
                [dist, dist]         => keep 1 trace at distance dist
                [distmin, distmax]   => keep distances in the range
                [d1, d2, d3, ... dN] => keep selected discrete distance
                default is None, keep all distances
        idrange = same as drange but given as index values
                [i, i]
                [imin, imax]
                [i1, i2, ..., iN]
        ddecim = 1 (default) or decimate

    return:
        an A1Section class instance containing header and data
    TODO implementer ddecim
    """
    import zmq
    import numpy as np
    import array
    from ._a1das_exception import ZMQSocketError, DataTypeError
    from sys import platform


    hd = a1.socket_header
    dhd = a1.data_header.copy()
    if hd is None:
        raise ZMQSocketError("socket not opened; open first before reading")

    # define block to read
    if block is None:
        nblock = 1
    elif isinstance(block, int):
        nblock=block
    elif len(block) < 1:
        nblock = 1
    elif len(block) == 1:
        nblock = max(int(block[0]), 1)
    elif len(block) == 2:
        nblock = block[1]-block[0]

    # select distance index
    if drange is not None:
        idrange = dhd.index(drange)
        idrange[1] += 1



    # windows/linux settings
    if platform == 'linux' or platform == 'darwin':
        idata_type = 'i'
    else:
        idata_type = 'l'

    # send ACQ
    try:
        hd.socket.send(np.double(0))
    except:
        raise ZMQSocketError("cannot send ACQ to "+hd.socketname+" socket")

    # receive header and data for 1st block
    try:
        message1 = hd.socket.recv()
        message2 = hd.socket.recv()
        message3 = hd.socket.recv()
    except:
        raise ZMQSocketError("cannot receive data from " + hd.socketname + " socket")

    #
    # read first block
    #

    DATA1 = array.array(idata_type, message1[0:4])
    count = DATA1[0]
    DATA1 = array.array('d', message1[4::])
    new_time_stamp = DATA1[0]
    DATA2 = array.array('d', message2[48:64])
    ntime = DATA2[1]-DATA2[0]
    #print('first block', count, new_time_stamp, hd.time_stamp)
    if count==4:
        type='raw'
    else:
        if ntime>0:
            type='strainrate'
        else:
            type='energyband'

    if (type!='strainrate'):
        raise DataTypeError('read socket only support strain rate')

    # reference block size for continuity checking
    block_time_size = dhd['ntime'] * 2
    block_space_size = dhd['nspace']

    #
    # make sure block was updated
    #
    old_time_stamp = hd.time_stamp
    while new_time_stamp == old_time_stamp:
        try:
            hd.socket.send(np.double(0))
            message1 = hd.socket.recv()
            message2 = hd.socket.recv()
            message3 = hd.socket.recv()
            DATA1 = array.array('d', message1[4::])
            new_time_stamp = DATA1[0]
            #print('time stamp',new_time_stamp,old_time_stamp)
        except:
            raise ZMQSocketError("cannot receive data from " + hd.socketname + " socket")
    #print('lu time stamp',new_time_stamp,old_time_stamp)
    a1.socket_header.time_stamp = new_time_stamp
    DATA2 = array.array(idata_type, message2[48:64])
    nspace = DATA2[1]-DATA2[0] + 1
    ntime = DATA2[3]-DATA2[2] + 1
    if ntime != block_time_size or nspace != block_space_size:
        print('ntime=',ntime,'block_time_size=',block_time_size,'nspace=',nspace,'bss=',block_space_size)
        raise ZMQSocketError("inconsistent block size with respect to previous read")

    #
    # build ndarray data block
    #
    DATA3 = array.array('f', message3)
    StrainRate = np.reshape(DATA3, (ntime, nspace))
    used_ntime = int(ntime/2)
    bs4 = int(ntime/4)
    if idrange is None:
        mspace = nspace
    elif len(idrange) == 2:
        mspace = len(range(idrange[0],idrange[1]))
    else:
        mspace = len(idrange)
    data = np.empty((used_ntime*nblock, mspace))

    if idrange is None:
        data[0:used_ntime,:]  = StrainRate[bs4:3*bs4, :]
    elif len(idrange) == 2:
        data[0:used_ntime, :] = StrainRate[bs4:3*bs4, idrange[0]:idrange[1]]
    elif len(drange) > 2:
        data[0:used_ntime, :] = StrainRate[bs4:3*bs4, idrange]

    #
    # Loop reading more blocks
    #
    for i in range(1, nblock):
        # read updated block
        while new_time_stamp == a1.socket_header.time_stamp:
            try:
                hd.socket.send(np.double(0))
                message1 = hd.socket.recv()
                message2 = hd.socket.recv()
                message3 = hd.socket.recv()
                DATA1 = array.array('d', message1[4::])
                new_time_stamp = DATA1[0]
            except:
                raise ZMQSocketError("cannot receive data from " + hd.socketname + " socket")

        a1.socket_header.time_stamp = new_time_stamp
        DATA2 = array.array(idata_type, message2[48:64])
        nspace = DATA2[1]-DATA2[0] + 1
        ntime = DATA2[3]-DATA2[2] + 1
        if ntime != block_time_size or nspace != block_space_size:
            print('ntime=',ntime,'block_time_size=',block_time_size,'nspace=',nspace,'bss=',block_space_size)
            raise ZMQSocketError("inconsistent block size with respect to previous read")

        #
        # build ndarray data block
        #
        DATA3 = array.array('f', message3)
        StrainRate = np.reshape(DATA3, (ntime, nspace))
        if idrange is None:
            data[i*used_ntime : (i+1)*used_ntime, :] = StrainRate[bs4:3*bs4, :]
        elif len(idrange) == 2:
            data[i*used_ntime  : (i+1)*used_ntime, idrange[0]:idrange[1]] = StrainRate[bs4:3*bs4,                                                       idrange[0]:idrange[1]]
        elif len(drange) > 2:
            data[i*used_ntime : (i+1)*used_ntime, idrange] = StrainRate[bs4:3*bs4, idrange]

    dhd.set_item(ntime=nblock * used_ntime)
    dhd.set_item(nspace=mspace)
    time = dhd['dt']*np.arange(bs4, 2*bs4*nblock + bs4) + new_time_stamp
    dhd.set_item(time = time)
    dhd.set_item(otime=time[0])

    return dhd, data


def close_das_socket(hd):
    '''

    '''
    hd.socket.close()

