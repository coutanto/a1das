#  Module to perform basic IO operations on A1 DAS data from file and socket stream,
#
# O. Coutant, ISTerre, UGA, 2020
#
# Content:
#   Definition of classes A1File, A1Section,...
#   (see details of content below)
#
# class A1File:
# ------------
#    members = file_header   (class _A1FileHeader)
#              socket_header (class _A1SocketHeader)
#              data_header   (dict)
#
#    methods:
#                __init__()
#                __str__() (print(A1File)
#
#               read()                    : read data with various selection options and return an A1Section instance
#               set_otime_from_filename() : set origin time from filename
#               close()                   : close file/socket stream
#    functions:
#               open()                    : open and return a file/socket description
#               tcp_address()             : format a socket address for ZMQ given IP address and port
#
# class A1Section:
# ---------------
#
#    members =  data_header   (dict)
#               data          (numpy 2D ndarray)
#    methods:
#               set_otime_from_filename() : set origin time from filename
#               obspy_stream()            : convert to obspy
#               time()                    : return time vector
#               dist()                    : return distance vector
#               plot()                    : simple plot vector
#               rplot()                   : simple plot ratser

#
# Content
#
# class _A1SocketHeader
# class _A1FileHeader
# "class" _A1DataHeader (actually not a class but just a dictionnary)
#
# class A1File
# class A1Section
# function open()
#
__version__ = "a1das Version 2.0.0"

__doc__="Generic IO functions to open/read/plot  febus/reducted/socket das files"
#
# ================================= _A1SocketHeader ================================
#
class _A1SocketHeader:
    def __init__(self,  socket=None,  socket_name=None, time_stamp=0., running_time=0.):
        self.socket = socket
        self.socketname = socket_name
        self.time_stamp = time_stamp
        self.running_time = running_time

    def __str__(self):
        from io import StringIO
        s = StringIO()
        s.write('\n<_A1SocketHeader>:')
        s.write('\n-----------')
        s.write('\nsocketname: '+self.socketname)
        s.write('\nrunning_time:'+str(self.running_time))
        return s.getvalue()

#
# ================================= _A1FileHeader ================================
#
class _A1FileHeader:

    def __init__(self,  freq_res, block_info, srate_node, chan_node, type, fd, filename):
        self.freq_res = freq_res    # frequency at which data are written on file
        self.srate_node = srate_node # hdf5 handle
        self.chan_node = chan_node   # hdf5 handle
        self.block_info = block_info # block info dict
        self.type = type                # 'febus' or 'reducted'
        self.fd = fd                    # file descriptor
        self.filename = filename        # as said

    def __str__(self):
        from io import StringIO
        s = StringIO()
        s.write('\n<_A1FileHeader>:')
        s.write('\n-----------------')
        s.write('\nfile type: '+self.type)
        s.write('\nfilename: ' + self.filename)
        if self.block_info is not None:
            for key, value in self.block_info.items():
                s.write('\n'+key+'='+str(value))
        return s.getvalue()


#
# ================================= class _A1DataHeader ================================
#
class _A1DataHeader:
    def __init__(self, gauge_length=None, sampling_res=None, prf=None, derivation_time=None,
                    data_type=None, data_axis=None, pos=None, dt=None, ntime=None, otime=None, dx=None, nspace=None,
                    ospace=None,dist=None, time=None, **kwargs):
        """
        fill the _A1DataHeader._header{} dictionnary with values passed as args

        gauge_length h   # gauge length
        sampling_res     # original spatial resolution
        derivation_time  # time step for time derivation
        prf              # pulse rate frequency
        pos              # 3D array of geographical coordinates
        data_type        # 'raw', 'strain', 'strain-rate'
        data_axis        # 'time_x_space' or 'space_x_time'
        dt               # time step
        ntime            # number of time samples
        time             # time vector
        otime            # origin time
        dx               # space step along fiber
        dist             # distance vector along fiber
        nspace           # number of space samples
        ospace           # origin position along fiber

        Note: why data_header is not a simple dictionnary?
        because we want to keep a predefined set of metadata, even if not defined (or None),
        and do this via the __init__() method that documents it
        """
        self._header = {
            'gauge_length': gauge_length,  # gauge length
            'sampling_res': sampling_res,  # original spatial resolution
            'derivation_time': derivation_time,  # time step for time derivation
            'prf': prf,  # pulse rate frequency
            'pos': pos,  # 3D array of geographical coordinates
            'data_type': data_type,  # 'raw', 'strain', 'strain-rate'
            'data_axis': data_axis,  # 'time_x_space' or 'space_x_time'
            'dt': dt,  # time step
            'ntime': ntime,  # number of time samples
            'otime': otime,  # origin time in Posix format
            'dx': dx,  # space step along fiber
            'nspace': nspace,  # number of space samples
            'ospace': ospace,  # origin position along fiber
            'time': time,
            'dist': dist
        }
        for key, value in kwargs.items():
            self._header[key] = value

    def __getitem__(self, item):
        """
        overload [] operator for the _A1DataHeader class
        """
        try:
            return self._header[item]
        except:
            return None

    def copy(self):
        """
        # Description
        Create and copy a data_header with a 'real' copy of the dictionnary
        #Description return:
        a new instance of _A1DataHeader filled with the value of the calling instance
        """
        dhdout = _A1DataHeader()
        dhdout._header = self._header.copy()
        return dhdout

    def __str__(self):

        from io import StringIO
        from numpy import ndarray
        from datetime import datetime, timezone

        unit = {'gauge_length': ' [m]',
                'sampling_res': ' [cm]',
                'derivation_time': ' [msec]',
                'dt': ' [sec]',
                'dx': ' [m]'
                }

        s = StringIO()
        s.write('\n<_A1DataHeader>:')
        s.write('\n-----------------')
        for key, value in self._header.items():
            if value is None:
                continue
            if isinstance(value, ndarray):
                continue
            try:
                su = unit[key]
                s.write('\n' + str(key) + '= ' + str(value) + su)
            except:
                s.write('\n' + str(key) + '= ' + str(value))

            if key == 'otime':
                s.write('  <=> date =' + str(datetime.fromtimestamp(value, tz=timezone.utc)))

        if self._header['time'] is not None:
            time = self._header['time']
            s.write('\n\nTotal time duration :' + str(time[-1] - time[0]) + ' [sec]')
        if self._header['dist'] is not None:
            dist = self._header['dist']
            s.write('\nTotal distance :' + str(dist[-1] - dist[0]) + ' [m]')

        return s.getvalue()

    def __repr__(self):
        return self.__str__()

    def set_item(self, **kwargs):
        for key, value in kwargs.items():
            self._header[key]=value

    def keys(self):
        return self._header.keys()
    #
    # ====================================   INDEX()  =======================================
    #
    def index(self, drange=None, trange=None):
        """
        ## Description
        Return the index or list of indices that correspond(s) to the list of distance(resp. time) range
        in the 'distance' vector (resp. time vector)

        ## Input
        !!!!! Only one of trange or drange can be given

        drange = [unique_dist]; [dist_min, dist_max]; [d1, d2, ... dN]  (default = None)

        trange = [unique_time]; [time_min, time_max]; [t1, t2, ... tN]  (default = None)

        ## Return
        A list of indices that matches the given range in the A1Section.data_header['dist'] or
        A1Section.data_header['time'] vectors

        ## Usage example
        >>> import a1das
        >>> f=open('filename')
        >>> a=f.read()
        >>> dhd = a.data_header
        >>> dlist1 = dhd.index(drange=[50.])          # single index at distance 50 m
        >>> dlist2 = dhd.index(drange=[50., 60.])     # 2 indices at 50 and 60m
        >>> dlist3 = dhd.index(drange=[[dist] for dist in np.arange(50,60.,1.)]) # all indices for distance 50< ..<60 every 1m
        """
        from numpy import nanargmin, abs, ndarray
        from ._a1das_exception import DataTypeError, WrongValueError

        if drange is not None and trange is not None:
            raise WrongValueError('cannot define trange AND drange, select one of two')

        if drange is not None:
            dtrange=drange
            vec = self['dist']
        else:
            dtrange=trange
            vec = self['time']

        if isinstance(dtrange, ndarray):
            raise DataTypeError('range must be a list or tuple')

        if dtrange is None:
            return [0, -1]

        d1 = nanargmin(abs(vec - dtrange[0]))
        indx = [d1]
        if len(dtrange) == 1:
            indx.append(d1+1)
        else:
            for i in range(1,len(dtrange)):
                indx.append(nanargmin(abs(vec-dtrange[i])))
        return indx


#
# ========================================= CLASS A1File ============================
#
class A1File:
    """
    ## Description
    A1File class contains file/socket descriptors to read and convert DAS data from Febus format,
    or H5 reducted file format or data socket stream (ZMQ protocol).

    ## Class content

        A1File.file_header:    file header    (_A1FileHeader class instance)
        A1File.socket_header:  socket header  (_A1SocketHeader class instance)
        A1File.data_header:    metadata       (_A1DataHeader class instance)

    The metadata list is obtained by typing the name of the A1File instance (see example below)

    ## Reading/setting data_header metadata
    Read directly as a dictionnary from a1_file (ex: f['dt']) or through data_header (ex: f.data_header['dt'])

    Set by the set_item method() ex: f.set_item(dt=0.1) or through data_header ex: f.data_header.set_item(dt=0.1)

    ## Usage example
        >>> import a1das
        >>> f = a1das.open("my_filename")
        >>> f # print metadata header values
        >>> GL = f['gauge_length']
        >>> f.set_item(gauge_length=10.)
    """
    def __init__(self, file_header=None, socket_header=None, data_header=None):
        # metadata
        self.file_header   = file_header   # _A1FileHeader class instance : io from file
        self.socket_header = socket_header # io from socket
        self.data_header   = data_header   # A1DasDataHeder class instance : some meta data

    def set_item(self, **kwargs):
        """
        ##Description
        Update/create a matadata entry in the data_header class dictionnary attached to the A1File

        ##Input
            **kwargs = arbitrary list of key_word = value
        ##Usage example
            >>> import a1das
            >>> f=a1das.open(filename)
            >>> f.set_item(dog_name='medor',dog_food='candies',dog_age=7)

        """
        self.data_header.set_item(kwargs)

    def __getitem__(self, item):
        """
        overload [] operator for the A1File class
        """
        try:
            return self.data_header._header[item]
        except:
            return None
    #
    # ========================================= READ() ============================
    #
    def read(self, block=None, trange=None, drange=None, ddecim=1, skip=True, verbose=0, **kwargs):
        """
        ## Description
        Read data from a Febus/reducted file or a socket directly connected to Febus DAS
        in memory and return a `A1Section` class instance.

        All data are read unless one of the three options is set: block, trange, drange.
        block is exclusive with trange,drange and return data aligned on block limits

        Data are read as float32 and are converted by default as float64 unless specified

        ## Input:
            block = (list), block=[start, end, step]. Read full block of data, from start to end by step.
                [0,0,?] read all blocks
                [] read one block
                default = None
            trange = (list), in sec, trange=[start_time, end_time], select time from start_time to end_time
                default = None
            drange = (list), in meter, drange=[start_position, end_position], select distance from start_position
                to end_position
                default = None
            ddecim = decimation along spatial dimension, (must be a divider of chunck/block size)

            skip = read all blocks/chunck of data (False) or one over 2 (True)
                blocks contains redundant data, skip=True is faster but may yield some errors on block limits

            verbose = verbosity level

            **kwargs= any optionnal supplementary arguments, see below

        ## Return:
        A `A1Section` class instance

        ## Usage example
            >>> import a1das
            >>> f = a1das.open("my_filename")
            >>> a1 = f.read(trange=(0.,1800.),drange=(100., 900.),skip=False)

        ## optionnal arguments
            float_type = 'float_32' or 'float_64' for febus file only
        """

        from ._core_febus_file import read_febus_file
        from ._core_reducted_file import read_reducted_file
        from ._core_socket import read_das_socket

        # read from a file (febus or reducted)
        if self.file_header is not None:
            if self.file_header.type == 'febus':
                data_header, data = read_febus_file(self, block, trange, drange, ddecim, skip, verbose=verbose,**kwargs)

            elif self.file_header.type == 'reducted':
                data_header, data = read_reducted_file(self, trange=trange, drange=drange, verbose=verbose,**kwargs)

        # read from a socket
        elif self.socket_header is not None:
            data_header, data = read_das_socket(self, block, drange=drange, ddecim=ddecim,**kwargs)

        else:
            data_header = data = None

        # create A1DSection class instance
        if self.file_header is not None:
            data_header.set_item(filename=self.file_header.filename)

        return A1Section(data_header, data)
    #
    # ========================================= CLOSE() ============================
    #
    def close(self):
        """
        Close a das file or socket stream

        """
        from ._core_febus_file import close_das_file
        from ._core_socket import close_das_socket

        # On lit depuis un fichier
        if self.file_header:
            close_das_file(self.file_header)
            self.file_header = None

        # On lit depuis un socket
        elif self.socket_header:
            close_das_socket(self.socket_header)
            self.socket_header = None


    #
    # ====================================    PRINT()  =======================================
    #
    def __str__(self, verbose=0):
        """
        Display the informations of a Febus DAS hdf5 file, called by print()

        verbose:  = if 1 (default) only acquisition information is printed; if >1, internal HDF5 file structure is printed
        """
        
        from io import StringIO
        s = StringIO()

        s.write('\n<A1File>:')
        s.write('\n-----------')
        if self.file_header is not None:
            s.write('\nfile_header:')
            # print file_header part
            s.write(self.file_header.__str__())

        if self.socket_header is not None:
            s.write('\nsocket_header:')
            # print socket_header part
            s.write(self.socket_header.__str__())

        # print data_header part
        s.write('\n\ndata_header:')
        s.write(self.data_header.__str__())
        return s.getvalue()

    def __repr__(self):
        """
        Display the informations of a Febus DAS hdf5 file, called by print()
        verbose:  = if 1 (default) only acquisition information is printed; if >1, internal HDF5 file structure is printed
        """
        return self.__str__()


    #
    # ====================================    _GET_TIME_BOUNDS()  =======================================
    #
    # Note O.C:
    #
    # Data are written by block (or chunk) @ Freq_res/1000 Hz (typically 1 Hz),
    # that contains several records
    # One record corresponds to all strain(rate) values along fiber recorded at a specific time
    # Each block contains "bts2" records of data (e.g 1sec) plus additional overlapping
    # The blocks overlap by half their length, (i.e. 1sec in the above example), they are thus 2sec long
    # One quarter (bts4) overlap with the previous block, one quarter (bts4) with the next one
    # The useful part of data, the <record> starts at 1/4 and ends at 3/4 of a block
    #
    #  [ - ~ ~ -]                    block 1
    #       [- ~ ~ -]                block 2
    #          [ - ~ ~ - ]           block 3
    #             .....
    #              [ - ~ ~ - ]       last block
    #
    #      - : unused duplicated time samples, there are bts4 time samples in '-'
    #      ~ : time samples read = <record>, there are bts2 time samples in '~ ~'
    #
    # block_time_size is the number of time samples in a block, bts2 is half of that, bts4 a quarter
    # dt is the time step between sample, i.e, records

    def _get_time_bounds(self, block=None, trange=None, skip=True, align_on_block=False, tdecim=None):
        """
        file_header._get_time_bounds(block, trange, skip, align_on_block, tdecim)

        Create block and time indices for reading according to the argument choice

        block    = [first_block, last_block] select data comprised between first_block, last_block (comprised),
                    default = None
        trange   = extract time starting at range[0] to range[1] or everything is trange=None and block = None
        skip     = skip redundant block when reading data (default = True)
        align_on_block = aligned data on start of block & end of block rather than exact time (if trange is set)
        tdecim   = time decimation factor

        return:  first_block,    = first block to be read in the file
                 last_block,     = last block +1 to be read in the file
                 step_block,     = step
                 time_vector_out = vector of absolute time value on output
                 time_indices    = [first, last] indices to be read in the 1st, last and other blocks

        """
        #
        from numpy import arange, fix

        nb_block_max = self.file_header.block_info['nb_block']  # total number of block in the file
        freq_res = self.file_header.freq_res  # frequency of block writing

        bts = self.file_header.block_info['time_size']
        bts2 = int(bts / 2)
        bts4 = int(bts / 4)
        dT = self.data_header['dt']
        origin_time = self.data_header['otime']

        # origin time of the file is the first sample read by febus (they skip one fourth of the first block)
        offset = bts4 * dT
        #offset_rule = offset
        # origin time of the file is the first sample of the file (this solution according to Gaëtan)
        offset_rule = 0.

        if skip:
            # !!!! if skip  = 2, one can read only an ODD number of blocks
            step_block = 2
            if nb_block_max % 2 == 0:
                nb_block_max = nb_block_max - 1  # The last EVEN block cannot be read
        else:
            # !!!! if skip  = 1, one can read any number of block
            step_block = 1
        #
        # define the read <for> loop indices: first_block, last_block
        #  for i, block in enumerate(range(first_block, last_block, step_block)):
        #
        if not trange:
            if block:
                first_block = block[0]
                last_block = min(block[1] + 1, nb_block_max)
                nb_block = last_block - first_block
            else:
                first_block = 0  # First block index for temporal acquisition
                nb_block = nb_block_max
                last_block = nb_block  # last block index in the reading loop, ie, never reached

            # for skip=2, read one block over 2, we scan an ODD number of block
            # start time is the beginning of the first block
            # !!!!! RECALL than nb_block IS NOT the number of blocks read, but the number of block scanned on the disk
            if skip:
                if nb_block % 2 == 0:
                    nb_block = max(nb_block + 1, nb_block_max)
                    last_block = first_block + nb_block - 1
                from_time = first_block * bts2 * dT
                # Total_time_size = number of time samples that we read in the time series
                total_time_size = (nb_block + 1) * bts2
                # total_time_size = nb_block * bts2 #ModifOC
                # time indices for 1st block, other blocks, last block
                time_indices = [[0, bts], [0, bts], [0, bts]]

            # for skip = 1, read every block, no restriction, but read only half of their size
            # start time is the first quarter of the first block because we don't want to manage
            # the 1st quarter of block (by laziness...)
            else:
                from_time = first_block * bts2 * dT + bts4 * dT
                total_time_size = nb_block * bts2
                # time indices for 1st block, other blocks, last block
                time_indices = [[bts4, 3 * bts4], [bts4, 3 * bts4], [bts4, 3 * bts4]]
        #
        # read only selected time range
        #
        else:
            from_time = trange[0]
            to_time = trange[1]

            if not skip:
                # check range that was given as argument
                min_time = bts4 * dT
                #min_time = 0.
                if from_time < min_time:
                    from_time = min_time
                max_time = ((nb_block_max * bts2 + bts4) * dT)
                if to_time > max_time:
                    to_time = (nb_block_max * bts2 + bts4) * dT
                    print('Warning: trange max is to high, set to:', str(to_time))

                # first_block = int(np.fix(from_time * freq_res))  # index of block containing the start time
                # last_block  = int(np.fix(to_time * freq_res)) + 1  # index of block containing the ending time + 1
                ##first_block = int(fix(from_time * freq_res)) - 1  # index of block containing the start time
                ##last_block = int(fix(to_time * freq_res))  # index of block containing the ending time + 1
                first_block = int(fix((from_time-offset) * freq_res))
                last_block  = int(fix((to_time-offset) * freq_res)) + 1
                if last_block > nb_block_max:
                    last_block = nb_block_max
                nb_block = last_block - first_block

                # round from_time to the start time of the first block that is read
                # time indices for 1st block, other blocks, last block
                # For future use: if we decide not to round to block limits

                if not align_on_block:
                    dt_1 = from_time - (first_block * bts2 + bts4) * dT
                    i_1 = int(dt_1 / dT)
                    dt_last = to_time - ((last_block - 1) * bts2 + bts4) * dT
                    i_last = int(dt_last / dT)
                    if nb_block >1:
                        time_indices = [[bts4 + i_1, 3 * bts4], [bts4, 3 * bts4], [bts4, i_last + bts4]]
                    else:
                        time_indices = [[bts4 + i_1, i_last + bts4], [bts4 + i_1, i_last + bts4], [bts4 + i_1, i_last + bts4]]
                    total_time_size = (nb_block * bts2) - i_1 - (bts2 - i_last)
                else:
                    time_indices = [[bts4, 3 * bts4], [bts4, 3 * bts4], [bts4, 3 * bts4]]
                    from_time = (first_block * bts2 + bts4) * dT
                    # Total_time_size = number of time samples in the time series
                    total_time_size = nb_block * bts2



            else:  # skip = True
                if from_time < 0:
                    from_time = 0
                max_time = ((nb_block_max + 1) * bts2 * dT)
                if to_time > max_time:
                    to_time = max_time

                # first_block = int(np.fix(from_time * freq_res))  # index of block containing the start time
                # last_block  = int(np.fix(to_time * freq_res)) + 1  # index of block containing the ending time
                ##first_block = int(fix(from_time * freq_res)) - 1 # index of block containing the start time
                ##last_block = int(fix(to_time * freq_res))  # index of block containing the ending time + 1
                first_block = int(fix((from_time) * freq_res / 2.))*2
                last_block  = int(fix((to_time) * freq_res / 2.))*2 + 1

                nb_block = last_block - first_block
                if nb_block % 2 == 0:  # Always read an ODD number of block
                    last_block += 1
                    nb_block += 1
                if last_block > nb_block_max:
                    last_block = nb_block_max

                # round from_time to the start time of the first block that is read
                # For future use: if we decide not to round to block limits
                if not align_on_block:
                    # starting time in the first block of the first point
                    dt_1 = from_time - first_block * bts2 * dT
                    i_1 = int(dt_1 / dT)
                    # ending time in the last block of the first point
                    dt_last = to_time - (last_block - 1) * bts2 * dT
                    i_last = int(dt_last / dT)
                    if nb_block > 1:
                        time_indices = [[i_1, bts], [0, bts], [0, i_last]]
                    else:
                        time_indices = [[i_1, i_last], [i_1, i_last], [i_1, i_last]]
                    total_time_size = (nb_block + 1) * bts2 - i_1 - (bts - i_last)
                else:
                    time_indices = [[0, bts], [0, bts], [0, bts]]
                    from_time = first_block * bts2 * dT
                    # Total_time_size = number of time samples in the time series
                    total_time_size = (nb_block + 1) * bts2

        # Create the time vector that contains samples datation
        # decimation can only be performed when reading entire blocks
        # The calling process must ensure that the block size is a multiplier of decimation factor
        # case no decimation
        if tdecim is None or tdecim == 1:
            time_vector_out = arange(0, total_time_size) * dT + from_time #+ origin_time
        # case no time range and time decimation
        elif not trange:
            time_vector_out = arange(0, total_time_size, tdecim) * dT + from_time #+ origin_time
        # case
        else:
            print('Time decimation is only valid when reading complete file or a block list')
            time_vector_out = arange(0, total_time_size) * dT + from_time #+ origin_time

        return first_block, last_block, step_block, time_vector_out, time_indices

    #
    # ====================================    _GET_SPACE_BOUND()  =======================================
    #
    def _get_space_bounds(self, drange=None, ddecim=None):
        '''
        A1File._get_space_bounds(drange, ddecim)

        :param drange: extract positions from from_to_position[0] to from_to_position[1], if not used set to None
        :param ddecim: space decimation step, if not used set to 1

        :return: distance_fiber_out, distance_indices, distance_fiber_in
                 distance_vector selected, array of distance indices in the input distance_vector
        '''
        from numpy import linspace, nanargmin

        ZI_start = self.data_header['ospace']
        dX = self.data_header['dx']
        distance_fiber_in = linspace(0, self.data_header['nspace'] - 1, self.data_header['nspace']) * dX + ZI_start
        #ZI_end = distance_fiber_in[-1]

        # extract everything
        if not drange:
            first_point = 0  # index for 1st point
            last_point = self.data_header['nspace'] - 1  # index for last point

        # extract a selected part
        else:
            # check range given as argument
            from_position = drange[0]
            to_position = drange[1]
            if from_position < distance_fiber_in[0]:
                from_position = distance_fiber_in[0]
            if to_position > distance_fiber_in[-1]:
                to_position = distance_fiber_in[-1]
                print('Warning: drange max is too high, set to ', str(to_position))

            # set indices
            first_point = int(nanargmin(abs(distance_fiber_in- from_position)))
            last_point  = int(nanargmin(abs(distance_fiber_in- to_position)))

        # distance_length = Total_space_size
        # distance_start = First_point
        # distance_end = distance_start + distance_length
        distance_fiber_out = distance_fiber_in[
                             first_point:last_point + 1:ddecim]  # add 1 because last index is excluded
        distance_indices = range(first_point, last_point + 1, ddecim)  # add 1 because last index is excluded

        return distance_fiber_out, distance_indices, distance_fiber_in

    #
    # ====================================   SET_OTIME_FROM_FILENAME()  =======================================
    #
    def set_otime_from_filename(self, offset=0., prefix=None):
        """
        ##Description
        Set data header origin time from the filename assuming Febus convention on the filename
        ex: SR_2021-08-26_14-32-39_UTC.h5.
        This call affect the origin time to file header and is propagated to all subsequent A1File.read() call
        This does not affect the values of the <time> vector which always refer to the origin (ie starting) time
        ##Input
            offset = a time offset to add/substract from the filename information
            prefix = a prefix ending by "_" in case the filename do not follow Febus convention (SR_, RAW_, S_, ...)
        """
        from datetime import timezone, datetime
        from ._a1das_exception import WrongValueError

        name = self.file_header.filename
        # expected to be of the form SR_2021-08-26_14-32-39_UTC.h5
        if prefix is None:
            start = name.find("_")
            ofs = 1
        else:
            start = name.find(prefix)
            ofs = len(prefix)

        if start == -1:
            raise WrongValueError('could not set origin time from filename, check filename and/or prefix')
            return
        end = name.find("_UTC")
        s = name[start+ofs:end]
        # convert from string
        try:
            d = datetime.strptime(s, "%Y-%m-%d_%H-%M-%S")
        except:
            raise ValueError('wrong date-time format, check file prefix')

        # convert to UTC
        dd = datetime(d.year, d.month, d.day, d.hour, d.minute, d.second, tzinfo=timezone.utc)
        # convert to Posix
        self.data_header.set_item(otime=dd.timestamp()+offset)

    def time(self):
        """
        ## Description
        Return the time vector defined in the A1File data_header, i.e. same as f['time']
        """
        return self.data_header['time']

    def dist(self):
        """
        ## Description
        Return the distance vector defined in the A1File data_header, i.e. same as f['dist']
        """
        return self.data_header['dist']

    #
    # ========================================= CLASS A1SECTION ============================
    #
class A1Section:
    """
    ## Description
    A1Section class contains das data and metadata read from file or TCP stream

    ## Class content

        A1Section.data:         2D float array (ndarray of shape [ntime x nspace] or [nspace x ntime])
        A1Section.data_header   metadata       (_A1DataHeader class instance)


    ## Reading/setting data_header metadata
    Read directly as a dictionnary from a1_section ex: a1['dt'] or through data_header ex: a1.data_header['dt']

    Set by the set_item method() ex: a1.set_item(dt=0.1) or through data_header ex: a1.data_header.set_item(dt=0.1)
    ##Usage example
        >>>import a1das
        >>>f=a1das.open(filename)   # open file
        >>>a1 = f.read()            # read section
        >>>a1                       # print header values
        >>>time_step=a1['dt']       # get time step
        >>>a1.set_item(dt=0.1)      # set time step

    """
    def __init__(self, dhd=None, data=None):
        self.data_header = dhd    # dictionnary of header values
        self.data = data          # numpy 2D ndarray.
                                  # If data_header['axis'] = 'time_x_space' first dimension is time
                                  # If data_header['axis'] = 'space_x_time' first dimension is space

    def time(self):
        """
        ## Description
        Return the time vector defined in the <A1Section> data_header, i.e. shortcut for a1.data_header['time']
        """
        return self.data_header['time']

    def dist(self):
        """
        ## Description
        Return the distance vector defined in the <A1Section> data_header, i.e. shortcut for a1.data_header['dist']
        """
        return self.data_header['dist']

    def otime(self):
        """
        ##Description
        Return the UTC origin (or start) time of the section
        in two formats: Posix time and string

        ##Return
        posix_otime, string_otime
        """
        from datetime import datetime, timezone
        time = self['time']   # could read otime too
        otime=time[0]

        return otime, datetime.fromtimestamp(otime, tz=timezone.utc).strftime('%Y:%M:%d-%H:%M:%S.%f')


    def __str__(self):
        from io import StringIO
        s = StringIO()

        s.write('\n<A1Section>:')
        s.write('\n--------------')
        s.write('\ndata: \nndarray data['+str(self.data_header['ntime'])+'x'+str(self.data_header['nspace'])+']\n')
        # print file_header part
        s.write('\ndata_header: ')
        s.write(self.data_header.__str__())

        return s.getvalue()

    def __repr__(self):
        return self.__str__()

    def set_item(self, **kwargs):
        """
        ##Description
        Update/create a matadata entry in the data_header class dictionnary attached to A1Section

        ##Input
            **kwargs = arbitrary list of key_word = value
        ##Usage example
            >>> import a1das
            >>> f=a1das.open(filename)
            >>> a=f.read()
            >>> a.set_item(cat_name='medor',cat_food='candies',cat_age=7)

        """
        self.data_header.set_item(kwargs)
        
    def __getitem__(self, item):
        """
        overload [] operator for the A1Section class
        """
        try:
            return self.data_header._header[item]
        except:
            return None
    #
    # ====================================   SET_OTIME_FROM_FILENAME()  =======================================
    #
    def set_otime_from_filename(self, offset=0., prefix=None):
        """
        ##Description
        Set data header origin time from the filename assuming Febus convention on the filename
        ex: SR_2021-08-26_14-32-39_UTC.h5.
        This call affect the origin time to file header and is propagated to all subsequent A1File.read() call
        This does not affect the values of the <time> vector which always refer to the origin (ie starting) time

        ##Input
            offset = a time offset to add/substract from the filename information
            prefix = a prefix ending by "_" in case the filename do not follow Febus convention (SR_, RAW_, S_, ...)
        """
        from datetime import timezone, datetime
        from ._a1das_exception import WrongValueError

        if 'filename' in self.data_header._header.keys():
            name = self.data_header['filename']
            # expected to be of the form SR_2021-08-26_14-32-39_UTC.h5
            if prefix is None:
                start = name.find("_")
                ofs = 1
            else:
                start = name.find(prefix)
                ofs = len(prefix)

            if start == -1:
                raise WrongValueError('could not set origin time from filename, check filename and/or prefix')
                return
            end = name.find("_UTC")
            s = name[start+ofs:end]
            # convert from string
            try:
                d = datetime.strptime(s, "%Y-%m-%d_%H-%M-%S")
            except:
                raise ValueError('wrong date-time format, check file prefix')

            # convert to UTC
            dd = datetime(d.year,d.month,d.day,d.hour,d.minute,d.second, tzinfo=timezone.utc)
            # convert to Posix
            self.data_header.set_item(otime=dd.timestamp()+offset)
        else:
            raise WrongValueError('<filename> is not defined in data header')

    #
    # ====================================   OBSPY_STREAM()  =======================================
    #
    def obspy_stream(self, drange=None, station_name=None):
        """
        ##Description
        Return an obspy stream from a DAS section with optional spatial range drange

        ##Input
            drange: = (dmin, dmax) list or tuple with minimal and maximal distance in meter
        ##Return
            An obspy stream
        """
        from obspy.core import Stream, Trace, UTCDateTime
        from ._a1das_exception import WrongValueError

        st = Stream()
        if drange is None:
            s0=0
            s1=self['nspace']
        else:
            s = self.index(drange)
            s0 = s[0]
            s1 = s[1]+1

        dist = self['dist']
        time = self['time']
        if station_name is None:
            station_name = 'DAS'
        for i in range(s0, s1):
            header = {'npts': self.data_header['ntime'],
                      'station': station_name,
                      'starttime': UTCDateTime(self.data_header['otime']+time[0]),
                      #'sampling_rate': 1./float64(self.data_header['dt']),
                      'location': str(int(dist[i]))+'m',
                      'channel': 'S',
                      'delta' : self.data_header['dt']
                      }
            if self.data_header['data_axis'] == 'time_x_space':
                st.append(Trace(data=self.data[:,i], header=header))
            elif self.data_header['data_axis'] == 'space_x_time':
                st.append(Trace(data=self.data[i,:], header=header))
            else:
                raise WrongValueError('<data_axis> header field is neither <time_x_space> nor <space_x_time>')
        return st

    # ====================================   INDEX()  =======================================
    #
    def index(self, drange=None, trange=None):
        """
        ## Description
        Return the index or list of indices that correspond(s) to the list of distance(resp. time) range
        in the 'distance' vector (resp. time vector)

        ## Input
        !!!!! Only one of trange or drange can be given

        drange = [unique_dist]; [dist_min, dist_max]; [d1, d2, ... dN]  (default = None)

        trange = [unique_time]; [time_min, time_max]; [t1, t2, ... tN]  (default = None)

        ## Return
        A list of indices that matches the given range in the A1Section.data_header['dist'] or
        A1Section.data_header['time'] vectors

        ## Usage example
        >>> import a1das
        >>> f=open('filename')
        >>> a=f.read()
        >>> dlist1 = a.index(drange=[50.])          # single index at distance 50 m
        >>> dlist1 = a.index(50.)                   # single index at distance 50 m
        >>> dlist2 = a.index(drange=[50., 60.])     # 2 indices at 50m and 60m
        >>> dlist3 = a.index(drange=[[dist] for dist in np.arange(50,60.,1.)]) # all indices for distance 50< ..<60 every 1m
        >>> tlist = a.index(trange=[10.,11.])       # time indices for time between 10. and 20s after 1st sample
        """
        from numpy import nanargmin, abs, ndarray
        from ._a1das_exception import DataTypeError, WrongValueError

        if drange is not None and trange is not None:
            raise WrongValueError('cannot define trange AND drange, select one of two')

        if drange is None and trange is None:
            raise WrongValueError('Please define one of drange(distance range) or trange (time range)')

        if drange is not None:
            dtrange=drange
            vec = self.data_header['dist']
        else:
            dtrange=trange
            vec = self.data_header['time']

        if isinstance(dtrange, ndarray):
            raise DataTypeError('range must be a list or tuple')

        if isinstance(dtrange,float):
            dtrange=[dtrange]

        if dtrange is None:
            return [0, -1]

        d1 = nanargmin(abs(vec - dtrange[0]))
        indx = [d1]
        if len(dtrange) == 1:
            indx.append(d1+1)
        else:
            for i in range(1,len(dtrange)):
                indx.append(nanargmin(abs(vec-dtrange[i])))
        return indx
    #
    # ====================================   PLOT()  =======================================
    #
    def plot(self, fig=None, clip=100, splot=(1, 1, 1), title='', max=100, amax=None, variable_area=False, drange=None, trange=None, redraw=None):
        """
        ##Description
        Produce a vector plot of the das section, optionnaly with variable area
        ##Input
            fig= figure handle (default None)
            clip= (int) clip amplitude at clip% of the max amplitude (default 100%)
            splot= (tupple) plot as a subplot (default subplot(1,1,1))
            title= (str) figure title (default None)
            max= (int) maximum number of traces on the plot (default 100)
            amax= maximum amplitude for clipping if not using clip percentage (default None)
            variable_area= (bool) plot with variable_are (half wiggle is filled in black) (default False)
            drange= (tupple) distance range [dmin, dmax]
            trange= (tupple) time range [tmin, tmax]
            redraw = {lines} redraw lines of a precedent figure  or None for new plot, requires fig argument

        ##Return
            figure_handle, curve_list

        ##Usage example
            >>>import a1das
            >>>f=a1das.open(filename)
            >>>a1=f.read()
            >>>a1.plot(clip=50,variable_area=True)
        """
        from a1das._plot import plot
        from ._a1das_exception import ReadDataError

        if self.data is None:
            raise ReadDataError('Hum Hum, read data before plotting ...')

        fig, redraw = plot(self,fig, clip, splot, title, max, amax, variable_area, drange, trange, redraw)

        return fig, redraw
    #
    # ====================================   RPLOT()  =======================================
    #
    def rplot(self, fig=None, cmap='RdBu', vrange=None, splot=(1, 1, 1), title='', drange=None, trange=None):
        """
        ##Description
        Produce a raster plot of the das section
        ##Input
            fig= figure handle (default None)
            cmap= colormap (default=RdBu)
            vrange= trace amplitude range (default None)
            splot= plot as a subplot (default subplot(1,1,1))
            title= figure title (default None)
            drange= distance range
            trange= time range

        ## Return
            figure_handle, axes_handle

        ##Usage example
            >>>import a1das
            >>>f=a1das.open(filename)
            >>>a1=f.read()
            >>>a1.rplot()
        """
        from a1das._plot import rplot
        from ._a1das_exception import ReadDataError

        if self.data is None:
            raise ReadDataError('Hum Hum, read data before plotting ...')

        fig, ax = rplot(self,fig, cmap, vrange, splot, title, drange, trange)

        return fig, ax

    def save(self,filename):
        """
        ## Description
        a.save(filename)

        Save a `A1Section` that has been read, extracted or converted from stream or files into the reducted file format.

        ## Input
        filename= a name for the hdf5 file that contains the reducted file

        ## Usage example
            >>> import a1das
            >>> f = a1das.open('my_filename', format='febus')
            >>> a = f.read(trange=(tmin,tmax), drange=(dmin,dmax), skip=False) # read and extract a subportion of the original file
            >>> a.save('new_section') # 'a' can be read later using the same a1das.open() and f.read() functions
            >>> f.close()

        """
        from ._core_reducted_file import _save_reducted_file
        _save_reducted_file(self,filename)

    def concat(self, b):
        """
        ## Description
            a.concat(b): Concatenate section b behind section a along the fast axis
            if a and b are ordered as [space x time], fast axis is time and on obtain a new section [ space x 2*time]
            if a and b are ordered as [time x space], not implemented

        ## Return
            a new section
        """
        import numpy as np
        from ._a1das_exception import DataTypeError

        if self['data_axis'] != 'time_x_space' or b['data_axis'] != 'time_x_space':
            raise DataTypeError('section must be ordered with data as time x space')

        time_a = self['time']
        time_b = b['time']

        nspace_a = self['nspace']
        nspace_b = b['nspace']
        if nspace_a != nspace_b:
            raise DataTypeError('concat: different spatial spacing')

        dta = self['dt']
        dtb = b['dt']
        if dta != dtb:
            raise DataTypeError('concat failed, different sampling rate')


        # Create new section
        data_header = self.data_header.copy()
        data_header.set_item(ntime=self['ntime']+b['ntime'])
        data_header.set_item(time=np.concatenate((time_a,time_b)))
        newdata = np.vstack((self.data, b.data))

        return A1Section(data_header, newdata)

#
# ========================================= OPEN()  ============================
#
def open(filename, format=None):
    """
    ## Description
    f=open('path',format)  where format = 'febus', 'reducted', 'socket'.

    Open a Febus das file, a reducted Hdf5 file or a socket stream on a DAS Febus interrogator and return an `A1File` instance

    ## Input
    filename= a filename or a socket address in the form "tcp://ip_address:port"

    format = one of 'febus', 'reducted', 'socket'. If not supplied, the function try to determine it automatically

    The socket address can be formatted using the a1das.tcp_address utility

    ## Return
    a `A1File` class instance

    ## Usage example
        >>> import a1das
        >>> f = a1das.open('my_filename', format='reducted')
        >>> f2 = a1das.open(a1das.tcp_address('128.0.0.1',6667), format='socket')

    """
    from ._core_febus_file import open_febus_file
    from ._core_reducted_file import open_reducted_file
    from ._core_socket import open_das_socket
    from ._a1das_exception import FileFormatError, WrongValueError

    if format == 'febus' or format is None:
        try:
            file_header, data_header = open_febus_file(filename)
            return A1File(file_header=file_header, data_header=data_header)
        except (KeyError, OSError, AttributeError):
            print('not in febus format')
            pass

    if format == 'reducted' or format is None:
        try:
            file_header, data_header = open_reducted_file(filename)
            return A1File(file_header=file_header, data_header=data_header)
        except (FileFormatError, OSError):
            print('not in reducted format')
            pass

    if format == 'socket' or format is None:
        try:
            socket_header, data_header = open_das_socket(filename)
            return A1File(socket_header=socket_header, data_header =data_header)
        except:
            raise FileFormatError('Cannot open <Febus>, <reducted> or <socket> file format>')

    # if we are here, some problem occured: wrong format argument for instance
    if format is not None:
        raise WrongValueError(" could not open file<"+filename+">, wrong format argument?")

def tcp_address(address='127.0.0.1', port=6667):
    """
    Return a string of the form <tcp://address:port> from an IP address and a port number
    input:
        address = IP adress as a string, default to '127.0.0.1'
        port = IP port default to 6667
    """
    return "tcp://%s:%d" % (address, port)
#
# ========================================= private stuff ===========================
#
def __visit_item__(name, node):
    """
    called by f.visititems()
    """
    import h5py

    if isinstance(node, h5py.Dataset):
        # node is a dataset
        print('DATASET : ', node.name)
        print('   Size :', node.size)
        print('   MaxSize : ', node.maxshape)
        print('   dtype : ', node.dtype)
        print('   chunkSize : ', node.chunks)

    else:
        # node is a group
        # list all attributes keys (it's a dict.)
        print('GROUP: ', node.name)
        print('Attributes: ')
        for key in node.attrs.keys():
            print('   ', key, node.attrs[key])
    #
    # display structure of hdf5 file on standard output
    #


def __display__(f):
    f.visititems(__visit_item__)


# Needed for compatibility with 2.7
# f.visit() can't call directly print in 2.7, but can call myprint()
def __myprint__(x):
    print(x)



