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
#              data_header   (class _A1DataHeader)
#
#    methods:
#                __init__()
#                __str__() (print(A1File)
#
#               read()                    : read data with various selection options and return an A1Section instance
#               close()                   : close file/socket stream
#    functions:
#               open()                    : open and return a file/socket description
#               tcp_address()             : format a socket address for ZMQ given IP address and port
#
# class A1Section:
# ---------------
#
#    members =  data_header   (class _A1DataHeader)
#               data          (numpy 2D ndarray)
#    methods:
#               obspy_stream()            : convert to obspy
#               time()                    : return time vector
#               dist()                    : return distance vector
#               sync()                    : sama as SAC, set time[0] to 0 and shift origin time
#               plot()                    : simple plot vector
#               rplot()                   : simple plot ratser
#
__version__ = "a1das Version 2.0.0"

__doc__="Generic IO functions to open/read/plot  febus/reducted/socket das files"

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

    def is_transposed(self):
        """
        ## Description
        Only relevant for 'reducted' format. Return True if the data are transposed from the original Febus format.

        Transposed =  the data array is stored as data[nspace, ntime],

        Not transposed = the data array is stored as data[ntime, nspace],

        ## Return
        False for Febus original files, True/False for reducted files according to f['axis?'] header value
        """
        if self['axis1'] == 'space':
            return False
        else:
            return True
    #
    # ========================================= READ() ============================
    #
    def read(self, block=None, trange=None, drange=None, skip=False, verbose=0, float_type='float64',
             tdecim=1, ddecim =1):
        """
        ## Description
        Read data from a Febus/reducted file or a socket directly connected to Febus DAS
        in memory and return a `A1Section` class instance.

        All data are read unless one of the three options is set: block, trange, drange.
        block is exclusive with trange,drange and return data aligned on block limits

        Data are read as float32 and are converted by default as float64 unless specified

        ## Input all format:
            trange = (list/tuple) time range in sec from 1st sample, [start, end]/[value] or
                     (range) index range  (default = None, read all)
            drange = (list/tuple) distance range in meter from 1st sample [start, end]/[value] or
                     (range) index range (default = None, read all)
            skip = (bool) read all blocks/chunck of data (False) or one over 2 (True) (default False)
                blocks contains redundant data, skip=True is faster but may yield some errors on block limits
            verbose = verbosity level
            float_type = 'float_32' or 'float_64' float type for data in memory, default is float64

        ## Input febus file
            block = (list), block=[start, end, step]. Read full block of data, from start to end by step.
                [0,0,?] read all blocks
                [] read one block
                default = None
            ddecim = (int) decimation along spatial dimension, (default 1, must be a divider of chunck/block size)

        ## Input reducted file and transposed data
            tdecim = (int) time decimation (default 1)

        ## Return:
        A `A1Section` class instance

        ## Usage example
            >>> import a1das
            >>> f = a1das.open("my_filename", format='febus')
            >>> a1 = f.read(trange=(0.,1800.),drange=(100., 900.),skip=False) #select by time and distance
            >>> a2 = f.read(trange=range[0,100]) #select by indices

        """

        from ._core_febus_file import read_febus_file
        from ._core_reducted_file import read_reducted_file
        from ._core_socket import read_das_socket

        #
        # read from a febus file file
        #
        if self.file_header is not None:
            if self.file_header.type == 'febus':
                data_header, data = read_febus_file(self, block = block, drange=drange, trange=trange, ddecim=ddecim,
                                                    skip = skip, verbose=verbose, float_type=float_type)
        #
        # read from a reducted file
        #
            elif self.file_header.type == 'reducted':
                data_header, data = read_reducted_file(self, drange=drange, trange=trange, verbose=verbose,
                                                       float_type=float_type, tdecim=tdecim)
        #
        # read from a socket
        #
        elif self.socket_header is not None:
            data_header, data = read_das_socket(self, block = block, drange=drange)

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
        ## Description
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
        a1file._get_time_bounds(block, trange, skip, align_on_block, tdecim)

        Create block and time indices for reading according to the argument choice

        block    = [first_block, last_block] select data comprised between first_block, last_block (comprised),
                    default = None
        trange   = [start, end] sec, extract time starting at trange[0] to trange[1] or everything is trange=None and block = None
                  start & end are counted with respect to file relative start-time (0=1st sample)
        skip     = skip redundant block when reading data (default = True)
        align_on_block = align data on start of block & end of block rather than exact time (if trange is set)
        tdecim   = time decimation factor

        return:  first_block,    = first block to be read in the file
                 last_block,     = last block +1 to be read in the file
                 step_block,     = step
                 time_vector_out = vector of absolute time value on output
                 time_indices    = [first, last] indices to be read in the 1st, last and other blocks

        """
        #
        from numpy import arange, fix
        from ._a1das_exception import WrongValueError

        nb_block_max = self.file_header.block_info['nb_block']  # total number of block in the file
        freq_res = self.file_header.freq_res  # frequency of block writing

        bts = self.file_header.block_info['block_time_size']
        bts2 = int(bts / 2)
        bts4 = int(bts / 4)
        dt = self.data_header['dt']
        origin_time = self.data_header['otime']

        # origin time of the file is the first sample read by febus (they skip one fourth of the first block)
        offset = bts4 * dt
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
                from_time = first_block * bts2 * dt
                # Total_time_size = number of time samples that we read in the time series
                total_time_size = (nb_block + 1) * bts2
                # total_time_size = nb_block * bts2 #ModifOC
                # time indices for 1st block, other blocks, last block
                time_indices = [[0, bts], [0, bts], [0, bts]]

            # for skip = 1, read every block, no restriction, but read only half of their size
            # start time is the first quarter of the first block because we don't want to manage
            # the 1st quarter of the block
            else:
                from_time = first_block * bts2 * dt + bts4 * dt
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
                min_time = bts4 * dt
                max_time = ((nb_block_max * bts2 + bts4) * dt)
                if from_time > max_time or to_time <min_time:
                    raise WrongValueError('time range requested is beyond time bounds')
                if from_time < min_time:
                    from_time = min_time
                    print('Warning: trange min is to small, set to:', str(from_time))
                if to_time > max_time:
                    to_time = (nb_block_max * bts2 + bts4) * dt
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
                    dt_1 = from_time - (first_block * bts2 + bts4) * dt
                    i_1 = int(dt_1 / dt)
                    dt_last = to_time - ((last_block - 1) * bts2 + bts4) * dt
                    i_last = int(dt_last / dt)
                    if nb_block >1:
                        time_indices = [[bts4 + i_1, 3 * bts4], [bts4, 3 * bts4], [bts4, i_last + bts4]]
                    else:
                        time_indices = [[bts4 + i_1, i_last + bts4], [bts4 + i_1, i_last + bts4], [bts4 + i_1, i_last + bts4]]
                    total_time_size = (nb_block * bts2) - i_1 - (bts2 - i_last)
                else:
                    time_indices = [[bts4, 3 * bts4], [bts4, 3 * bts4], [bts4, 3 * bts4]]
                    from_time = (first_block * bts2 + bts4) * dt
                    # Total_time_size = number of time samples in the time series
                    total_time_size = nb_block * bts2



            else:  # skip = True
                min_time = 0.
                max_time = ((nb_block_max + 1) * bts2 * dt)
                if from_time > max_time or to_time <min_time:
                    raise WrongValueError('time range requested is beyond time bounds')
                if from_time < min_time:
                    from_time = 0
                    print('Warning: trange min is to small, set to:', str(from_time))
                if to_time > max_time:
                    to_time = max_time
                    print('Warning: trange max is to high, set to:', str(to_time))
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
                    dt_1 = from_time - first_block * bts2 * dt
                    i_1 = int(dt_1 / dt)
                    # ending time in the last block of the first point
                    dt_last = to_time - (last_block - 1) * bts2 * dt
                    i_last = int(dt_last / dt) + 1 #modifOC
                    if nb_block > 1:
                        time_indices = [[i_1, bts], [0, bts], [0, i_last]]
                    else:
                        time_indices = [[i_1, i_last], [i_1, i_last], [i_1, i_last]]
                    total_time_size = (nb_block + 1) * bts2 - i_1 - (bts - i_last)
                else:
                    time_indices = [[0, bts], [0, bts], [0, bts]]
                    from_time = first_block * bts2 * dt
                    # Total_time_size = number of time samples in the time series
                    total_time_size = (nb_block + 1) * bts2

        # Create the time vector that contains samples datation
        # decimation can only be performed when reading entire blocks
        # The calling process must ensure that the block size is a multiplier of decimation factor
        # case no decimation
        if tdecim is None or tdecim == 1:
            time_vector_out = arange(0, total_time_size) * dt + from_time
        # case no time range and time decimation

        elif not trange:
            time_vector_out = arange(0, total_time_size, tdecim) * dt + from_time
        # case
        else:
            #print('Time decimation is only valid when reading complete file or a block list')
            time_vector_out = arange(0, total_time_size, tdecim) * dt + from_time

        return first_block, last_block, step_block, time_vector_out, time_indices

    #
    # ====================================    _GET_SPACE_BOUND()  =======================================
    #
    def _get_space_bounds(self, drange=None, ddecim=None):
        '''
        A1File._get_space_bounds(drange, ddecim)

        :param drange: extract positions from from_to_position[0] to from_to_position[1], if not used set to None
        :param ddecim: space decimation step, if not used set to 1

        :return:
        distance_fiber_out: vector of distances requested by drange and decime arguments
        distance_indices:   range of indices to build distance_fiber_out
        distance_fiber_in:  original vector of distances

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
                          or
                                4 x 2D int8 array (ndarray of shape |4 x ntime x nspace]) for febus RAW data
        A1Section.data_header   header       (_A1DataHeader class instance)


    ## Reading/setting data_header key+value
    The list of header (key,value) store in a A1section named 'a' is obtained by typing: print(a)

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

    def __str__(self):
        """
        invoked by print()
        """
        from io import StringIO
        s = StringIO()

        s.write('\n<A1Section>:')
        s.write('\n--------------')
        s.write('\ndata: \nndarray data[' + str(self.data.shape[0]) + 'x' + str(self.data.shape[1]) + ']')

        # print file_header part
        s.write('\ndata_header: ')
        s.write(self.data_header.__str__())

        return s.getvalue()

    def __repr__(self):
        return self.__str__()

    def __getitem__(self, *items):
        """
        overload [] operator for the A1Section class
        a['field'] returns the header field value
        a[n1:n2,m1:m2] returns the data slices
        """
        # request for header value
        from numpy import nanargmin
        if len(items) == 1 and isinstance(items[0],str):
            return self.data_header._header[items[0]]

        slices = items[0]
        # request for data slice

        if isinstance(slices[0],slice) and isinstance(slices[1],slice):
            islice1, islice2 = self.data_header.dim2index(slices[0],slices[1])
            return self.data[islice1, islice2]

    #
    #
    #
    def subsection(self,trange=(None,None),drange=(None,None)):
        """
        ## Description
        return a A1Section which is a subsection of the original section

        ## Input
        trange = (list) [start, stop] time range in sec starting from origin time. (default = (None, None), keep all)
        drange = (list) [start, stop]] distance range in m starting from origin. (default = (None, None), keep all)
        """
        tslice=slice(trange[0],trange[1],None)
        dslice=slice(drange[0],drange[1],None)
        if self['axis1'] == 'time':
            itslice, idslice = self.data_header.dim2index(tslice, dslice)
        else:
            idslice, itslice = self.data_header.dim2index(dslice, tslice)

        data_header = self.data_header.copy()

        data_header.set_item(dist=self['dist'][idslice])
        data_header.set_item(nspace=len(self['dist'][idslice]))
        data_header.set_item(time=self['time'][itslice])
        data_header.set_item(ntime=len(self['time'][itslice]))

        if self['axis1'] == 'time':
            newdata = self.data[itslice, idslice].copy()
        else:
            newdata = self.data[idslice, itslice].copy()

        return A1Section(dhd=data_header, data=newdata)
    #
    # ====================================   time()  =======================================
    #
    def time(self):
        """
        ## Description
        Return the time vector defined in the <A1Section> data_header, i.e. shortcut for a1.data_header['time'] or a1['time']
        """
        return self.data_header['time']
    #
    # ====================================   dist()  =======================================
    #
    def dist(self):
        """
        ## Description
        Return the distance vector defined in the <A1Section> data_header, i.e. shortcut for a1.data_header['dist'] or a1['dist']
        """
        return self.data_header['dist']
    #
    # ====================================   otime()  =======================================
    #
    def otime(self):
        """
        ##Description
        Return the UTC origin (or start) time of the section
        in two formats: Posix time and string

        ##Return
        posix_otime, string_otime
        """
        from datetime import datetime, timezone
        otime = self['otime']

        return otime, datetime.fromtimestamp(otime, tz=timezone.utc).strftime('%Y:%m:%d-%H:%M:%S.%f')

    #
    # ====================================   set_item()  =======================================
    #
    def set_item(self, **kwargs):
        """
        ##Description
        Update/create a header entry in the data_header class dictionnary attached to A1Section

        ##Input
            **kwargs = arbitrary list of key_word = value
        ##Usage example
            >>> import a1das
            >>> f=a1das.open(filename)
            >>> a=f.read()
            >>> a.set_item(cat_name='medor',cat_food='candies',cat_age=7)

        """
        self.data_header.set_item(kwargs)

    #
    # ====================================   is_transposed()  =======================================
    #
    def is_transposed(self):
        """
        ## Description
        Only relevant for <reducted> format, return True if the data are transposed in file from the original Febus format.
        Transposed =  the data array is stored as data[nspace, ntime],
        Not transposed = the data array is stored as data[ntime, nspace],

        ## Return
        False for Febus original files, True/False for reducted files
        """
        # test on first axis <time> or <space>
        if self['axis1'] == 'time':
            return False
        else:
            return True
    #
    # ====================================   Synchronize()  =======================================
    #
    def synchronize(self):
        """
        ## Description
        This function is similar to the synchronize command in SAC.
        Set the origin time (field 'otime' in data header) to the time of the first sample
        in the section. After this operation, the value of time[0] is equal to 0.
        """
        otime = self['otime']
        time = self['time']
        otime += time[0]
        time -= time[0]
        self.set_item(otime=otime)
        self.set_item(time=time)

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
            if self.data_header['axis1'] == 'time':
                st.append(Trace(data=self.data[:,i], header=header))
            elif self.data_header['axis1'] == 'space':
                st.append(Trace(data=self.data[i,:], header=header))

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

        drange = (list or tuple) [unique_dist]; [dist_min, dist_max]; [dist1, dist2, ... distN] (default = None, take all)

        trange = (list or tuple) [unique_time]; [time_min, time_max]; [t1, t2, ... tN]  (default = None, take all)

        ## Return
        A list of indices that match the given range in the A1Section.data_header['dist'] or
        A1Section.data_header['time'] vectors

        ## Usage example
        >>> import a1das
        >>> f=open('filename')
        >>> a=f.read()
        >>> dhd = a.data_header

        >>> dlist1 = dhd.index(drange=50.)            # get index for distance 50.m
        >>> dlist1 = dhd.index(drange=[50.])          # get index at distance 50 m
        >>> dlist2 = dhd.index(drange=[50., 60.])     # get 2 indices at 50 and 60m
        >>> dlist3 = dhd.index(drange=[[dist] for dist in np.arange(50,60.,1.)]) # get all indices for distance 50< ..<60 every 1m
        >>> dlist4 = dhd.index(drange=(50.,60.))        # get all indices between 50. and 60.m
        >>> tlist = dhd.index(trange=[10., 20.])      # 2 indices for time 10s and 20 sec after 1st sample
        """
        return self.data_header.index(drange=drange, trange=trange)

    #
    # ====================================   PLOT()  =======================================
    #
    def plot(self, fig=None, clip=100, splot=(1, 1, 1), title='', max=100, amax=None, by_trace=False, variable_area=False, drange=None, trange=None, redraw=None):
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
            by_trace = (bool) normalize by plot (False) or by trace (True) (default=False)
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
        fig, redraw = plot(self,fig, clip, splot, title, max, amax, by_trace, variable_area, drange, trange, redraw)

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
    #
    # ====================================   save()  =======================================
    #
    def save(self,filename):
        """
        ## Description
        a.save(filename)

        Save a `A1Section` in a file using the reducted file format.

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
    #
    # ====================================   concat()  =======================================
    #
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

        if self['axis1'] !=  b['axis1']:
            raise DataTypeError('section must be ordered with data as time x space')

        dta = self['dt']
        dtb = b['dt']
        if dta != dtb:
            raise DataTypeError('concat failed, different sampling rate')

        time_a = self['time']
        time_b = b['time'] + time_a[-1] + dta
        print()

        nspace_a = self['nspace']
        nspace_b = b['nspace']
        if nspace_a != nspace_b:
            raise DataTypeError('concat: different spatial spacing')



        # Create new section
        data_header = self.data_header.copy()
        data_header.set_item(ntime=self['ntime']+b['ntime'])
        data_header.set_item(time=np.concatenate((time_a,time_b)))
        newdata = np.hstack((self.data, b.data))

        return A1Section(data_header, newdata)

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
            if format is None:
                print('open Febus format')
            return A1File(file_header=file_header, data_header=data_header)
        except (KeyError, OSError, AttributeError):
            pass

    if format == 'reducted' or format is None:
        try:
            file_header, data_header = open_reducted_file(filename)
            if format is None:
                print('open reducted format')
            return A1File(file_header=file_header, data_header=data_header)
        except (FileFormatError, OSError):
            pass

    if format == 'socket' or format is None:
        try:
            socket_header, data_header = open_das_socket(filename)
            if format is None:
                print('open ZMQ socket format')
            return A1File(socket_header=socket_header, data_header =data_header)
        except:
            raise FileFormatError('Cannot open <Febus>, <reducted> or <socket> file format>')

    # if we are here, some problem occured: wrong format argument for instance
    if format is not None:
        raise WrongValueError(" could not open file<"+filename+">, wrong format argument?")

#
# ========================================= tcp_address()  ============================
#
def tcp_address(address='127.0.0.1', port=6667):
    """
    Return a string of the form <tcp://address:port> from an IP address and a port number
    input:
        address = IP adress as a string, default to '127.0.0.1'
        port = IP port default to 6667
    """
    return "tcp://%s:%d" % (address, port)
#
# ========================================= set_true_time_policy()  ============================
#
def set_true_time_policy(policy=True):
    """
    ##Description
    set_true_time_policy()

    Define if the true time policy is applied (True) or not (False)

    When the true time policy is used, the timing of a sample is defined from the start time
    of the block it was contained in.

    When the true time policy is not used, one assume a constant time step for the time serie

    Default is true_time_policy=True
    """
    from ._core_febus_file import _set_true_time_policy
    _set_true_time_policy(policy)
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



