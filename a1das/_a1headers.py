#
#
# Content: definitions of the header classes: data header, file header and socket header
#
#
#          Although the data header contains just a dictionary, it is defined as a class
#          to have a constructor that defines precisely a (limited) number of keys
#


#
# Definition of required header keys and their meaning
#
_required_header_keys = {
    'gauge_length':      'gauge length (meter)',
    'sampling_res':      'original spatial resolution (cm), usually between 20 and 60 cm',
    'prf':               'laser Pulse Rate Frequency (Hz)',
    'data_type':         'data type = raw, strain, strain-rate, ...',
#    'data_axis':         'data axis phys. dim. : [time, space] (orig.) or [space, time] (transposed)',
    'axis1':             'dimension of data first axis: "time" or "space"',
    'axis2':             'dimension of data second axis: "time" or "space"',
    'dt':                'time step (sec) for time axis',
    'ntime':             'number of time samples',
    'otime':             'origin time. Absolute time of a sample is : time + otime',
    'dx':                'spatial step (m) for distance axis',
    'nspace':            'number of space samples',
    'ospace':            'origin position along fiber. The absolute distance is given by dist + ospace',
    'dist':              'relative distance vector w/r to origin position',
    'time':              'relative time vector w/r to origin time'
}
# Definition of other header keys that can be filled on the fly
_other_header_keys = {
    'derivation_time':       'time interval (sec) used for time derivation',
    'time_derivation_order': 'finite difference order for time derivation',
    'space_derivation_order': 'finite difference order for space derivation',
}

#
# ================================= class _A1DataHeader ================================
#
class _A1DataHeader:
    def __init__(self, hdr=None):
        """
        Fill the _A1DataHeader._header{} dictionnary with dictionary passed as arg
        It must contains the required information as definied in _required_header_keys dict

        """
        from ._a1das_exception import MissingHeaderKeywordError

        self._header={}

        # 1) handle the transposition info that has changed from old format
        try:
            val = hdr['transposition']
            if val == 1:
                hdr['axis1'] = 'space'
                hdr['axis2'] = 'time'
            else:
                hdr['axis2'] = 'space'
                hdr['axis1'] = 'time'
            del hdr['transposition']
        except:
            pass

        # 2) check required (mandatory) header keys
        for keys in _required_header_keys:
            if not keys in hdr:
                raise MissingHeaderKeywordError('Missing header mandatory key: ' + keys)

        # 3) copy keys
        self._header = { key: val for key, val in hdr.items()}
        #for keys in hdr: TODO à virer si correct
        #    self._header[keys] = hdr[keys]

    def __getitem__(self, item):
        """
        overload [] operator for the _A1DataHeader class
        """
        #try:
        return self._header[item]
        #except:
            #return None

    def copy(self):
        """
        # Description
        Create and copy a data_header with a 'real' copy of the dictionnary
        #Description return:
        a new instance of _A1DataHeader filled with the value of the calling instance
        """
        #dhdout = _A1DataHeader()
        #dhdout._header = self._header.copy()
        dhdout = _A1DataHeader(hdr = self._header)
        return dhdout

    def __str__(self):
        """
        specify the print() output
        """
        from io import StringIO
        from numpy import ndarray
        from datetime import datetime, timezone

        s = StringIO()
        s.write('\n<_A1DataHeader>:')
        s.write('\n-----------------')
        for key, value in self._header.items():
            if value is None:
                continue
            if isinstance(value, ndarray):
                s.write('\n' + str(key) + ' = ndarray of shape :'+str(value.shape))
                continue
            try:
                #su = unit[key]
                su = _required_header_keys[key]
                s.write('\n' + str(key) + '= ' + str(value) + "\t| "+su)
            except:
                try:
                    su = _other_header_keys[key]
                    s.write('\n' + str(key) + '= ' + str(value) + "\t| "+su)
                except:
                    s.write('\n' + str(key) + '= ' + str(value))

            if key == 'otime':
                s.write('\notime date (UTC) = ' + str(datetime.fromtimestamp(value, tz=timezone.utc)))

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
        """
        method to define or change the value of a key in the header
        """
        for key, value in kwargs.items():
            self._header[key]=value

    def keys(self):
        """
        Return the list of keys of the header dict
        """
        return self._header.keys()

    def items(self):
        """
        Return the pairs {keys, value} for the header dict
        """
        return self._header.items()

    def dim2index(self,slice1,slice2):
        """
        return the index slice given a dimension slice.
        i.e. convert from dimension to index into data array
        """
        from numpy import nanargmin

        slices=(slice1,slice2)
        if self['axis1'] == 'time':
            it = 0
            id = 1
        else:
            it = 1
            id = 0
        # time
        time = self['time']
        tstart = tend = tstep = None
        if slices[it].start is not None:
            tstart = nanargmin(abs(time - slices[it].start))
        if slices[it].stop is not None:
            tend = nanargmin(abs(time - slices[it].stop))
        if slices[it].step is not None:
            tstep = int(slices[it].step / self['dt'])
        if tend == tstart:
            tend += 1
        # space
        dist = self['dist']
        dstart = dend = dstep = None
        if slices[id].start is not None:
            dstart = nanargmin(abs(dist - slices[id].start))
        if slices[id].stop is not None:
            dend = nanargmin(abs(dist - slices[id].stop))
        if slices[id].step is not None:
            dstep = int(slices[id].step / self['dx'])
        if dstart == dend:
            dend += 1
        #print(dstart, dend, dstep,slices[id])
        if self['axis1'] == 'time':
            return slice(tstart, tend, tstep), slice(dstart, dend, dstep)
        else:
            return slice(dstart, dend, dstep), slice(tstart, tend, tstep)
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

        >>> dlist1 = dhd.index(drange=50.)            # get index for distance 50.m
        >>> dlist1 = dhd.index(drange=[50.])          # get index at distance 50 m
        >>> dlist2 = dhd.index(drange=[50., 60.])     # get 2 indices at 50 and 60m
        >>> dlist3 = dhd.index(drange=[[dist] for dist in np.arange(50,60.,1.)]) # get all indices for distance 50< ..<60 every 1m
        >>> dlist4 = dhd.index(drange=(50.,60.))        # get all indices between 50. and 60.m
        >>> tlist = dhd.index(trange=[10., 20.])      # 2 indices for time 10s and 20 sec after 1st sample
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

        if isinstance(dtrange, float) or isinstance(dtrange, int):
            dtrange = [float(dtrange)]

        if isinstance(dtrange, ndarray):
            raise DataTypeError('t(d)range must be a list, a tuple or a range')

        if dtrange is None:
            return [0, -1]

        d1 = nanargmin(abs(vec - dtrange[0]))
        indx = [d1]
        if len(dtrange) == 1:
            indx.append(d1+1)

        elif isinstance(dtrange,list):
            for i in range(1,len(dtrange)):
                indx.append(nanargmin(abs(vec-dtrange[i])))

        elif isinstance(dtrange,tuple):
            d2 = nanargmin(abs(vec - dtrange[1]))
            for i in range(d1,d2+1):
                indx.append(i)
        return indx

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

    def __init__(self,  freq_res, block_info, node, type, fd, filename):
        self.freq_res = freq_res     # frequency at which data are written on file
        self.node = node             # hdf5 handle to data
        self.block_info = block_info # block info dict
        self.type = type                # 'febus' or 'reducted'
        self.fd = fd                    # file descriptor
        self.filename = filename        # as said

    def __str__(self):
        from io import StringIO
        from numpy import ndarray
        s = StringIO()
        s.write('\n<_A1FileHeader>:')
        s.write('\n-----------------')
        s.write('\nfile type: '+self.type)
        s.write('\nfilename: ' + self.filename)
        if self.block_info is not None:
            for key, value in self.block_info.items():
                if isinstance(value, ndarray):
                    s.write('\n' + str(key) + ' = ndarray of shape :' + str(value.shape))
                    continue
                s.write('\n'+key+'='+str(value))
        return s.getvalue()
