from a1das._a1das_exception import VirtualServerError
import array
import numpy as np
import struct

__doc__="Functions to emulate a DAS interrogator socket server. The DAS is replaced by a file"+ \
        "\n These functions allows to play file data, or any process applied on data, and view them like a movie"

def ZMQVirtualDasServer(file, host='127.0.0.1', port=6667, block=None):
    """
    ## Description
    Read Febus strain[rate] HDF5 file and send to a ZMQ socket stream

    ## Input
        file = hdf5 (.h5) file to read
        host = IP adress for the virtual server
        port = IP port for the virtual server

    ## Usage example
        >>> from a1das import vserver
        >>> vserver.ZMQVirtualDasServer('my_SR_file',host='127.0.0.1',port=6667)
    """
    import a1das
    import zmq
    import pickle
    from io import BytesIO
    from timeit import default_timer as timer

    # open ZMQ connection
    context = zmq.Context()
    socket = context.socket(zmq.REP)
    socket.bind("tcp://%s:%d" % (host, port))
    print('launching DAS virtual server on address: <',"tcp://%s:%d" % (host, port),'>')

    #
    # open file
    #
    f = a1das.open(file)
    # set origin time
    #f.set_otime_from_filename()
    #read block info
    nblock = f.file_header.block_info['nb_block']
    #read freq_res
    freq_res = f.file_header.freq_res
    delta_t = 1./freq_res
    #time = f.data_header['otime']
    block_times = f.file_header.block_info['block_times']
    time = block_times[0]


    # wait for an acq
    _wait_for_acq(socket, time)

    # Set and Send header(s)
    message2 = _set_and_send_header_to_ZMQ(f.data_header, f.file_header, time, socket)

    data = a1das._core_febus_file._read_febus_file_block(f, 0)
    time0 = timer()
    socket.send(data)

    if block is None:
        block = [0, nblock]

    for i in range(block[0], block[-1]):
        #
        # wait for acq
        #
        _wait_for_acq(socket, time)

        #
        # resend headers and data
        #
        time1 = timer()
        #temporisation
        #while (time1 - time0 < delta_t):
        #    time1 = timer()
        time0 = time1
        #time += delta_t
        time = block_times[i]
        message2 = _set_and_send_header_to_ZMQ(f.data_header, f.file_header, time, socket, message2=message2)

        data = a1das._core_febus_file._read_febus_file_block(f, block=i)
        socket.send(data)


    socket.close()
    f.close()


#
# ========================================= ZMQRAWVIRTUALDASSERVER()  ============================
#
def ZMQRawVirtualDasServer(filein,  GL, DT, order_time=2, order_space=2, trange=None, drange=None, tdecim=1, ddecim=1,
                verbose=0, host='127.0.0.1', port=6667):
    """
    ## Description
    Convert Febus RAW data HDF5 file to strain[rate] and send them to a ZMQ socket stream

    ## input:
        filein: hdf5 (.h5) file to read
        GL: Gauge length (in meter)
        DT: Derivation time (IN SECOND)
        order_time:  finite derivation order in time, no derivation if set to 0 (default 2)
        order_space: finite derivation order in space (default 2)
        trange:  time range in sec (start,end), (default = None, everything is read)
        drange:  distance range in meter (not km) (start,end), (default = None, everything is read)
        ddecim:  distance decimation factor (default=1)
        tdecim:  time decimation factor (default=1)
        verbose: be verbose (default=0, minimum message level)


    ## Usage example
        >>> from a1das import vserver
        >>> gauge_length=10.
        >>> time_derivation=0.005
        >>> vserver.ZMQRawVirtualDasServer('my_RAW_file',gauge_length, time_derivation, host='127.0.0.1',port=6667)
    """

    import h5py
    from .core import open
    from a1das import _raw2strPy
    from ._a1das_exception import WrongValueError
    from timeit import default_timer as timer
    import zmq

    kchunk = 1
    skip = False
    host = '127.0.0.1'
    port = 6667
    use_compression = False
    #
    #  --------------------------- open file for reading
    #
    a1 = open(filein, format='febus')

    # open ZMQ connection
    context = zmq.Context()
    socket = context.socket(zmq.REP)
    socket.bind("tcp://%s:%d" % (host, port))

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
        exit(0)

    #
    #   --------------------------- compute time bounds and time vector -----------------
    #

    # indices for chunk (block) of data
    # indices for first and last time record in the first and last block
    # Vector of times in range[from_time, to_time] with tdecim
    # (block, trange, skip, align_on_block, tdecim)
    first_block, last_block, step_block, \
    time_out, block_indices = a1._get_time_bounds(trange=trange, skip=skip)

    #
    #     ---------------------------  compute distance bounds and indices ----------------
    #
    dist_out, dist_ix, dist_in = a1._get_space_bounds(drange, ddecim)
    if drange is not None:
        ospace_drange = drange[0]
    else:
        ospace_drange = 0.

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

    ncblock = kchunk * tdecim


    #
    # create messages from header values
    #
    freq_res = a1.file_header.freq_res
    delta_t = 1./freq_res
    time = dhd['otime']

    #
    # test GL and DT parameters
    #
    if GL < dhd['sampling_res'] / 100:
        raise WrongValueError(
            'Gauge length is smaller than sampling resolutionv' + str(GL) + '<' + str(dhd['sampling_res']))
    if DT < (1. / dhd['prf']):
        raise WrongValueError('time derivation is smaller than 1/Pulse_Rate_Frequency')

    #
    # --------------------------- loop reading blocks ----------------------
    #

    buff_out = np.empty((time_chunk_size, output_space_size, 4), np.int8, 'C')
    strain = np.empty((time_chunk_size, output_space_size), np.float64, 'C')
    print('size of data is '+str(time_chunk_size)+' x '+str(output_space_size) )

    time0 = timer()
    from timeit import default_timer as timer
    last_block_read = list(range(first_block, last_block, step_block))[-1]
    message2 = None
    for i, block in enumerate(range(first_block, last_block, step_block)):

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


        #
        # We have filled the temporary buffer
        # extract time decimated part and convert to strain
        #
        if jchunk == ncblock - 1:

            phase1_r = np.squeeze(buff_out[0:block_time_size * ncblock:tdecim, :, 0]).astype('int32')
            phase1_i = np.squeeze(buff_out[0:block_time_size * ncblock:tdecim, :, 1]).astype('int32')
            phase2_r = np.squeeze(buff_out[0:block_time_size * ncblock:tdecim, :, 2]).astype('int32')
            phase2_i = np.squeeze(buff_out[0:block_time_size * ncblock:tdecim, :, 3]).astype('int32')
            strain, GL, DT = _raw2strPy.raw2strpy(strain, phase1_r, phase1_i, phase2_r, phase2_i, dist_out, time_out,
                                                 GL, DT, order_space, order_time, verbose)
            data = strain.astype('float32')
            time1 = timer()
            while (time1 - time0 < delta_t):
                time1 = timer()
            time0 = time1

            _wait_for_acq(socket,time)
            message2 = _set_and_send_header_to_ZMQ(dhd, hd, time, socket, message2=message2)
            socket.send(data)
            time += delta_t



    # end of file reached, write partially filled buffer
    if jchunk != ncblock - 1:
        phase1_r = np.squeeze(buff_out[0:(jchunk + 1) * block_time_size:tdecim, :, 0])
        phase1_i = np.squeeze(buff_out[0:(jchunk + 1) * block_time_size:tdecim, :, 1])
        phase2_r = np.squeeze(buff_out[0:(jchunk + 1) * block_time_size:tdecim, :, 2])
        phase2_i = np.squeeze(buff_out[0:(jchunk + 1) * block_time_size:tdecim, :, 3])
        strain, GL, DT = _raw2strPy.raw2strpy(strain, phase1_r, phase1_i, phase2_r, phase2_i, dist_out, time_out,
                                             GL, DT, order_space, order_time, verbose)
        data = strain.astype('float32')
        # write to socket
        _wait_for_acq(socket,time)
        message2 = _set_and_send_header_to_ZMQ(dhd, hd, time, socket, message2=message2)
        socket.send(data)

    a1.close()
    socket.close()

#
# ========================================= ZMQREDUCTEDVIRTUALDASSERVER()  ============================
#
def ZMQReductedVirtualServer(filein, prefix, trange=None, drange=None, tdecim=1, ddecim=1,
                verbose=0, host='127.0.0.1', port=6667, block_time=2.):
    """
    ##Description
    Read a reducted file and send data on a ZMQ socket server
    ##Input
    filein = (string) file name
    prefix = (string) prefix to be used to extract origin time from filename
    trange = [tmin, tmax] time range
    drange = [dmin, dmax] distance range
    tdecim = (int) time decimation
    ddecim = (int) space decimation
    verbose = (int) verbosity level
    host = (string) IP address of host
    port = (int) port
    block_time = length of time block sent
    """
    import a1das
    import zmq
    import pickle
    from io import BytesIO
    from timeit import default_timer as timer

    f = a1das.open(filein, format='reducted')

    #get header usefull values
    dt = f['dt']
    nspace = f['nspace']
    ntime = f['ntime']

    #set distance reading limits
    if drange is not None:
        dlist = f.index(drange=drange)
        space_range = range(dlist[0],dlist[1],ddecim)
    else:
        space_range = range(0,nspace,ddecim)

    # set time reading limits
    if trange is not None:
        tlist = f.index(trange=trange)
        time_range = range(tlist[0], tlist[1],tdecim)
    else:
        time_range=(0,ntime,tdecim)

    dt = dt*tdecim

    #set Febus-file-like parameters
    block_time_size = int(block_time / dt)
    shift_time_size = int(block_time_size/2)
    delta_t = shift_time_size*dt
    f.file_header.block_info={'block_time_size':block_time_size}

    # set origin time
    #f.set_otime_from_filename(prefix=prefix)
    time = f['otime']



    # open ZMQ connection
    context = zmq.Context()
    socket = context.socket(zmq.REP)
    socket.bind("tcp://%s:%d" % (host, port))
    print('launching DAS virtual server on address: <',"tcp://%s:%d" % (host, port),'>')

    # wait for an acq
    _wait_for_acq(socket, time)

    # Set and Send header(s)
    message2 = _set_and_send_header_to_ZMQ(f.data_header, f.file_header, time, socket)

    #
    # Loop reading data
    #
    if f.is_transposed():

        #get HDF5 file descriptor
        fd = f.file_header.fd

        buffer = np.empty((len(space_range),block_time_size),dtype=np.float32)
        # Loop on block time
        for i in range(0,ntime,shift_time_size):
            start = i
            end = i + block_time_size
            for j, jx in enumerate(space_range):
                buffer[j,:] = fd['/Traces/'+str(jx)][start:end]

            message3 = buffer.transpose().copy()
            socket.send(message3)

            time += delta_t
            _wait_for_acq(socket, time)
            message2 = _set_and_send_header_to_ZMQ(f.data_header, f.file_header, time, socket, message2=message2)
    else:
        print('not yet implemented')

    f.close()
    socket.close()

def _wait_for_acq(socket, time=0):
    """
    Wait for an ACQ message from the server
    """
    from datetime import datetime, timezone
    msg = socket.recv()
    acq = struct.unpack("@d", msg)[0]
    if acq != 0.:
        raise VirtualServerError('Wrong Acq')
    else:
        print('received ', acq,' @ time',time,datetime.fromtimestamp(time, tz=timezone.utc).strftime('%Y:%m:%d-%H:%M:%S.%f'))

def _set_and_send_header_to_ZMQ(dhd, fhd, time, socket, message2=None):
    """
    The protocol contains 3 messages, 2 for headers and one for data
    This routine create/update message 1, prepare message 2 and send message1 and message2
    """
    import zmq
    from io import BytesIO

    message1 = BytesIO()
    message1.write(np.int32(1))
    message1.write(np.float64(time))

    if message2 is None:
        message2 = BytesIO()
        message2.write(np.float64(dhd['dx']))
        message2.write(np.float64(dhd['dt']*1000.))
        message2.write(np.float64(dhd['ospace']))
        message2.write(np.float64(0.))
        message2.write(np.float64(0.))
        message2.write(np.float64(0.))
        message2.write(np.int32([0]))
        message2.write(np.int32(dhd['nspace']-1))
        message2.write(np.int32([0]))
        message2.write(np.int32(fhd.block_info['block_time_size']-1))
        message2.write(np.float64(0.)) # to fill 72 bytes
        message2.write(np.int32([0]))

    socket.send(message1.getvalue(),flags=zmq.SNDMORE)
    socket.send(message2.getvalue(),flags=zmq.SNDMORE)

    return message2