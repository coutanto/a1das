#
# Module to read HDF5 files created by Febus A1 DAS interrogator
# O. Coutant, ISTerre, UGA, 2020 from matlab code by Febus
#

#
# Content:
#            Class Section
#
#
#            Functions:
#         s=a1das.load (filename) : load class instance from disk
#         a1das.info(filename)
#
# 1.0.1: typos sur l'appel du constructeur dans read et subsection
#        conflit de nom sur l'attribut infos, changé en attribut log
# 1.1.1: ajout des fonctions reduction et reduction_mpi
#        reduction_mpi necessite mpi4py, uen version de l'api hdf5 compilé avec l'option //
#        et h5py construit avec cette version hdf5
# 1.1.2: reduction2_mpi, reduction3_mpi, reduction4_mpi differents algos possibles
# 1.1.3: suppression des version mpi, version avec filtrage openMP via le paquet sosfilter, lecture 1 bloc / 2
#
import numpy as np
import h5py
import pickle
import scipy.io
import matplotlib.pyplot as plot

__version__ = "a1das Version 1.1.3"


class DasFileHeader:
    '''
    Class describing the metadata read from das hdf5 file header

    gauge_length (in meter)
    block_time_size : data are written by chunks that are block_time_size long
    nb_block :  number of data chunks
    freq_res : a time chunk is written every 1/freq_res sec
    dist_info : information about fiber distance,  dictionnary with {'step','npts', 'start', 'end'}
    time_info : information about time,            dictionnary with {'step','npts', 'start', 'end'}
    srate_node : hfd5 node or handle to read data
    '''

    def __init__(self, gauge_length, block_time_size, nb_block, freq_res, dist_info, time_info, srate_node):
        self.gauge_length = gauge_length
        self.block_time_size = block_time_size
        self.nb_block = nb_block
        self.freq_res = freq_res
        self.dist_info = dist_info
        self.time_info = time_info
        self.srate_node = srate_node

    @classmethod
    def read(self, f, Zone=1):
        '''
        Read the das file header and fill in the DasFileHeader class
        :param f: hdf5 file pointer
        :return: DasFileHeader class instance
        '''
        #
        # Organisation:    / DAS / SOURCE / ZONE
        #

        # Initialisation
        # Récupération des données d'enetete, dans les attributs du groupe Source1
        #
        print('')
        print('> Initialization')
        DAS = list(f.keys())  # only 1 elemts in this dict.
        DASName = DAS[0]  # 1st key gives name of groupe
        DASNode = f['/' + DASName]
        SOURCES = list(DASNode.keys())
        SOURCEName = SOURCES[0]
        SOURCENode = f['/' + DASName + '/' + SOURCEName]
        SrcAttr = SOURCENode.attrs  # it's a dict.
        ZONEName = 'Zone' + str(Zone)
        ZONENode = f['/' + DASName + '/' + SOURCEName + '/' + ZONEName]
        ZneAttr = ZONENode.attrs

        # DAS_name = DASName
        # Frequency = SrcAttr['PulseRateFreq'] / 1000.
        freq_res = SrcAttr['FreqRes'] / 1000.  # A chunk(block) of data is written @ (Freq_res/1000) Hz
        # Acquisition_length = SrcAttr['FiberLength']
        # Sampling_Res = SrcAttr['SamplingRes']

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

        gauge_length = ZneAttr['GaugeLength']

        TIMEName = 'time'
        # TIMENode = f['/' + DASName + '/' + SOURCEName + '/' + TIMEName]

        #
        # block structure and length
        #
        SRATEName = '/' + DASName + '/' + SOURCEName + '/' + ZONEName + '/' + 'StrainRate'
        SRATENode = f[SRATEName]
        # chunkSize = SRATENode.chunks
        nb_block = SRATENode.shape[0]
        block_time_size = SRATENode.shape[1]
        block_space_size = SRATENode.shape[2]

        #
        # distance information
        #
        Extent = ZneAttr['Extent']
        # ZI_start = Origin[0] + Extent[0] * dx
        # ZI_end = Origin[0] + Extent[1] * dx
        dist_info = {'step': dx, 'npts': block_space_size, 'start': Origin[0] + Extent[0] * dx, 'end': Origin[0] +
                                                                                                       Extent[1] * dx}

        #
        # time information
        #
        origin_time = Origin[1]  # WARNING WARNING actually this is not set in the header, timing could be obtained
        # from it,e.g. Posix time of the recording start
        time_info = {'step': dt, 'npts': nb_block, 'start': origin_time,
                     'end': (nb_block * int(block_time_size / 2)) * dt + origin_time}

        return DasFileHeader(gauge_length, block_time_size, nb_block, freq_res, dist_info, time_info, SRATENode)


class DasSection:
    """
    Class describing a <seismic> section recorded by Febus A1 DAS and read from their HDF5 file format
    Attributes:

      t         = time vector
      d         = distance vector
      data      = 2D data matrix (ndarray), strain rate (or strain ?)
      transpose = False if fast axis is distance (default), True if fast axis is time
      data      = 2D data vector,
                  [Time, distance] the fast dimension is the distance if transpose=False
                  [distance, time] the fast dimension is the time if transpose=True
      pos[3,:]  = 3D positions vector (optional)
      gauge_length = ...
      type      = 'strain rate' or 'strain'
      log       = dictionary with various information, initial file and later processing

    """

    #
    # =========================================== constructor ================================
    #
    def __init__(self, t=None, d=None, data=None, gauge_length=None, transpose=None, type=None, pos=None, log=None):
        '''
        :param t:
        :param d:
        :param data:
        :param gauge_length:
        :param transpose:
        :param type:
        :param pos:
        '''
        self.t = t
        self.d = d
        self.data = data
        self.gauge_length = gauge_length
        self.transpose = transpose
        self.type = type
        self.pos = pos
        self.log = log

    #
    # ====================================    READ()  =======================================
    #
    # Note for beginners like me: force classmethod so we can call it without a class instance
    #
    @classmethod
    def read(self, filename, trange=None, drange=None, tdecim=1, ddecim=1, Zone=1, Source=1):
        '''
        read a DAS section from an hdf5 Febus-Optics format
        ---------------------------------------------------

        :param filename: hdf5 (.h5) file to read
        :param trange:  time range in sec (start,end), (default = [] everything is read)
        :param drange:  distance range in meter (not km) (start,end), (default = [] everything is read)
        :param ddecim:  distance decimation factor (default=1)
        :param tdecim:  time decimation factor (default=1)
        :param Zone:    not used yet, (default=1)
        :param Source:  not used yet, (default=1)

        return a DasSection object
        '''

        #
        #  --------------------------- Time Range option -----------
        #

        if trange:
            section_time = True
            from_time = trange[0]  # in sec
            to_time = trange[1]
        else:
            section_time = False
            from_time = None
            to_time = None

        #
        #  --------------------------- distance Range option -----------
        #

        if drange:
            section_space = True
            from_position = drange[0]
            to_position = drange[1]
        else:
            section_space = False
            from_position = None
            to_position = None

        #
        #  --------------------------- read header ----------------------
        #
        #  open file
        f = h5py.File(filename, 'r')

        hd = DasFileHeader.read(f)

        #
        #   --------------------------- compute time bounds and time vector -----------------
        #

        # indices for chunk (block) of data
        # indices for first and last time record in the first and last block
        # Vector of times in range[from_time, to_time] with tdecim
        first_block, last_block, step_block, \
        time_out, time_ix = \
            __get_time_bounds__(hd, from_time, to_time, tdecim, section_time)

        #
        #     ---------------------------  compute distance bounds and indices ----------------
        #

        dist_out, dist_ix = __get_space_bounds__(hd, from_position, to_position, ddecim, section_space)

        #
        # ---------------------------   print summary  -------------------
        #
        print(' ')
        print('> Data extraction from [', time_out[0], ' - ', time_out[-1], '] sec and from [', dist_out[0], ' - ',
              dist_out[-1], '] m')

        #
        # --------------------------- initialize data ndarray ---------------------------
        #
        total_time_size = len(time_out)
        total_space_size = len(dist_out)
        StrainRate = np.empty((total_time_size, total_space_size), np.float32, 'C')

        #
        # --------------------------- loop reading blocks ----------------------
        #
        time_out_offset = 0
        for block in range(first_block, last_block):
            print('    ' + str(block - first_block + 1) + '/' + str(last_block - first_block) + ' blocks', end='\r')

            # compute time indices to be read in current block
            #time_ix = __get_time_indices_in_block__(block, first_block, last_block, hd.block_time_size, first_record,
            #                                       last_record, tdecim)
            time_length = len(time_ix)

            # copy current block data into section ndarray
            StrainRate[time_out_offset:time_out_offset + time_length, :] = hd.srate_node[block,
                                                                  time_ix[0]:time_ix[-1]+1:tdecim,
                                                                  dist_ix[0]:dist_ix[-1]+1:ddecim]
            time_out_offset += time_length

        f.close()

        #
        # --------------------------- set information in info dictionary ---------------------------
        #
        log = {'filename': filename, 'from_time': from_time, 'to_time': to_time, 'from_Position': from_position,
               'to_Position': to_position, 'tdecim': tdecim, 'ddecim': ddecim}

        #
        # ---------------------------  call constructor and return a DasSection instance --------------------------
        #
        return DasSection(t=time_out, d=dist_out, data=np.nan_to_num(StrainRate),
                          gauge_length=hd.gauge_length, transpose=False, type='strain rate',
                          pos=np.empty([3, dist_out.size], np.float64), log=log)

    #
    # ====================================    SUBSECTION()  =======================================
    #
    def subsection(self, trange=None, drange=None, tdecim=1, ddecim=1, copy=True):
        '''
        Extract a sub-section from a full section and/or decimate along one or two axis
        -------------------------------------------------------------------------------
        !!!! WARNING !!!!  if a subsection is created not in copy mode (copy=False)
        then all modifications are also applied to the original data.
        E.G. If you decimate in time and apply a filter, then original data will be screwed up
        RULE: 1) if you create a subsection for later processing, use true copy (default, copy=True)
              2) if you create a subsection without changing data (e.g. plotting), force copy=False

        trange = range for time extraction [min, max], (default=None)
        drange = range for distance extraction [min, max], (default=None)
        tdecim = time decimation (default=1)
        ddecim = distance decimation (default=1)
        copy   = create a new section with new memory allocation (copy=True) or
                 just follow a link toward original data
        '''

        if self.transpose == True:
            print('subsection method not implemented for transposed section')
            return None

        # convert to local
        if not drange:
            drange = [self.d[0], self.d[-1]]
        else:
            drange = list(drange)  # make sure drange is a list, not a tupple
        if not trange:
            trange = [self.t[0], self.t[-1]]
        else:
            trange = list(trange)

        # check range
        if drange[0] < self.d[0]:
            drange[0] = self.d[0]
        if drange[1] > self.d[-1]:
            drange[1] = self.d[-1]
        if trange[0] < self.t[0]:
            trange[0] = self.t[0]
        if trange[1] > self.t[-1]:
            trange[1] = self.t[-1]

        # steps
        dd = self.d[1] - self.d[0]
        dt = self.t[1] - self.t[0]
        # start indices
        id = int(np.fix((drange[0] - self.d[0]) / dd))
        it = int(np.fix((trange[0] - self.t[0]) / dt))
        # end indices
        jd = int(np.fix((drange[1] - self.d[0]) / dd))
        jt = int(np.fix((trange[1] - self.t[0]) / dt))

        print('subsection defined between [', self.d[id], self.d[jd], '] m and [', self.t[it], self.t[jt], '] sec')

        if copy:
            self.log.update({'subsection_trange': trange, 'subsection_drange': drange,
                             'subsection_tdecim': tdecim, 'subsection_ddecim': ddecim})
        if copy:
            return DasSection(t=self.t[it:jt + 1:tdecim].copy(), d=self.d[id:jd + 1:ddecim].copy(),
                              data=self.data[it:jt + 1:tdecim, id:jd + 1:ddecim].copy(), gauge_length=self.gauge_length,
                              transpose=False, type=self.type, pos=None, log=self.log)
        else:
            return DasSection(t=self.t[it:jt + 1:tdecim], d=self.d[id:jd + 1:ddecim],
                              data=self.data[it:jt + 1:tdecim, id:jd + 1:ddecim],
                              gauge_length=self.gauge_length, transpose=False, type=self.type, pos=None, log=self.log)

    #
    # ====================================    DUMP()  =======================================
    #
    def dump(self, filename, format='python'):
        '''
        Dump the instance of the "Section" class into a python binary file ".obj" if type = 'python' (default)
        or into a matlab ".mat" file if type='matlab'
        The python file will be read by load()

        format = format (default='python')
        '''
        if format == 'python':
            if ".obj" in filename:
                pickle.dump(self, open(filename, 'wb'))
            else:
                pickle.dump(self, open(filename + '.obj', 'wb'))

        elif format == 'matlab':
            if not ".mat" in filename:
                filename = filename + '.mat'
            scipy.io.savemat(filename, mdict={'D': self.d, 'T': self.t, 'Transp': self.transpose, 'DATA': self.data,
                                              'file': self.log['filename']},
                             do_compression=True)

    # =================================================================================
    # ====================================  plotting utilities ========================
    # =================================================================================
    @classmethod
    def fig(self, size=(18, 9)):
        '''
        Initialiaze a plot
        :param size: (width, heigth) (default=(18,9)
        :return: fig handle
        '''
        return plot.figure(figsize=size)

    @classmethod
    def show(self):
        '''
        show plot on screen
        '''
        plot.ion()
        plot.show()

    @classmethod
    def plotpause(self):
        '''
        When ran from ide or interpretor, figure is/are erased asap script ends
        pause to watch the plot
        '''
        plot.pause(0.001)
        input("Press [enter] to continue.")

    #
    # ====================================    RPLOT()  =======================================
    #
    def rplot(self, fig=None, cmap='RdBu', vrange=None, splot=(1, 1, 1), title=''):
        '''
        plot a DAS section as a raster image
        ------------------------------------
        fig = a fig handle necessary to plot several plot on the same figure (default=None, a new figure is created)
        cmap = colormap (default='RdBu')
        vrange = (vmin,vmax) range for colormap (default=[] set to max value)
        splot = subplot position (default=(1,1,1))
        '''

        if not fig:
            fig = plot.figure(figsize=(18, 9))

        if not vrange:
            vmax = np.max(np.absolute(self.data), axis=None)
            vmin = -vmax
        else:
            vmin = vrange[0]
            vmax = vrange[1]

        ax = plot.subplot(splot[0], splot[1], splot[2])
        extent = [self.d[0], self.d[-1], self.t[-1], self.t[0]]
        pos = ax.imshow(self.data, cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto', extent=extent)

        ax.set_xlabel('distance meters')
        ax.set_ylabel('time sec')
        ax.set_title('section DAS ' + title)
        fig.colorbar(pos)

    #
    # ====================================    PLOT()  =======================================
    #
    def plot(self, fig=None, clip=100, splot=(1, 1, 1), title='', max=100, amax=None, variable_area=False):
        '''
        plot a DAS section as traces
        ----------------------------
        fig   = a fig handle necessary to plot several plot on the same figure (default=None, a new figure is created)
        clip  = % of max amplitude that is clipped (deault=100) (clip and amax are exclusice)
        amax  = maximum amplitude that is clipped (default=None)
        splot = subplot position (default=(1,1,1))
        max=  = max number of trace per plot along distance dimension (default=100)
        variable_area = a kind of variable area plot, only positive wiggles are plotted (default=False)
        '''

        if not fig:
            plot.figure(figsize=(18, 9))

        from_distance = self.d[0]
        to_distance = self.d[-1]

        # Do not plot more than MAX traces
        # adjust decimation rate on distance
        ndist = len(self.d)  # number of distance
        if ndist < max:
            ddecim = 1
        else:
            ddecim = int(np.fix(ndist / max))

        idrange = list(range(0, ndist, ddecim))
        ntraces = len(idrange)
        max_value = np.max(np.absolute(self.data), axis=None)
        clip = clip / 100
        if not amax:
            clipped_value = max_value * clip
        else:
            clipped_value = amax
        width = 2 * clipped_value * ntraces
        offset = np.arange(0, ntraces + 1) * 2 * clipped_value
        gain = (to_distance - from_distance) / width
        if gain == 0:
            gain = 1

        # compute the width, trace offset, and gain

        ax = plot.subplot(splot[0], splot[1], splot[2])
        for i in range(0, len(idrange)):
            v = self.data[:, idrange[i]].copy()
            v[v > clipped_value] = clipped_value
            if not variable_area:
                v[v < -clipped_value] = -clipped_value
            else:
                v[v < 0] = 0
            v = v + offset[i]
            v *= gain
            plot.plot(self.t, v)

        ax.set_xlabel('time sec')
        ax.set_title('section DAS ' + title)

    #
    # ====================================  end of plotting utilities ========================
    #

    #
    # ====================================    SETPOSITION()  =======================================
    #
    def setPosition(self, filename, type='xyz'):
        '''
        Read position of the das line from a file and set pos[3,...] ndarray, possible types are
        'xyz', 'gpx', 'kml'
        The two extreme positions are supposed to be the fiber extremities,
        the positions are then interpolated along the fiber shape using cubic spline
        :param filename:
        :param type:
        :return:
        '''
        if type == 'xyz':  # read positions from a text file
            print('not yet implemented')
        elif type == 'gpx':
            print('not yet implemented')
        elif type == 'kml':
            print('not yet implemented')


#
# ==================================== end of class methods ============================================
#

def __visit_item__(name, node):
    '''
    called by f.visititems()
    '''
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

#
# ====================================    GET_TIME_BOUNDS()  =======================================
#
def __get_time_bounds__(hd, from_time, to_time, tdecim, section_time, skip=False):
    '''
    get_time_bounds(hd, from_time, to_time, tdecim, section_time, skip=False)

    :param : hd =  a class DasFileHeader instance
    :param from_time = extract time starting at <from_time>
    :param to_time   = until                   <to_time>
    :param tdecim:   = decimation factor
    :param section_time = False if we keep all time window, True otherwise
    :param skip = skip redundant block when reading data (default = False)
    :return: first_block,   = first block to be read in the file
             last_block,    = last block to be read in the file
             time_vector,   = vector of absolute time value in the file
             time_indices   = sample indices to be read in the blocks
                              We read the same for all blocks, rounding
                              to the closest second the requested time
    '''

    dT = hd.time_info['step']               # sample step in time
    origin_time = hd.time_info['start']     # absolute origin time
    nb_block = hd.nb_block                  # total number of block in the file
    freq_res = hd.freq_res                  # block writing frequency

    bts2 = int(hd.block_time_size / 2)
    bts4 = int(hd.block_time_size / 4)

    if skip:
        freq_res = 2 * freq_res
        block_time_size = hd.block_time_size
        step_block = 2
    else:
        block_time_size = bts2
        step_block = 1
    #
    # read all data
    #
    if not section_time:
        first_block = 0  # First block index for temporal acquisition
        last_block = nb_block  # last block index ....
        # Total_time_size = number of time samples in the time series
        total_time_size = nb_block * bts2 # this value is correct regardless of the choice for skip
        from_time = origin_time

    #
    # read only selected data
    #
    else:
        from_time -= origin_time  # WARNING WARNING assuming that Origin[1] contains starting time
        to_time -= origin_time

        # check range that was given as argument
        if from_time < 0:
            from_time = 0
        if to_time > ((nb_block * bts2 ) * dT):
            to_time = (nb_block * bts2 ) * dT
            print('Warning: trange max is to high, set to:', str(to_time))

        first_block = int(np.fix(from_time * freq_res))  # index of block containing the start time
        last_block = int(np.fix(to_time * freq_res))  # index of block containing the ending time
        if last_block > nb_block: # This can arise when skip=True
            last_block = nb_block

        # Total_time_size = number of time samples in the time series
        total_time_size = len(range(first_block, last_block+1, step_block)) * block_time_size

    # Create the time vector that contains sample datation
    time_vector = np.arange(0, total_time_size) * dT + from_time  # ModifOC

    # if decimation requested, we need to recompute
    # number of samples and time vector
    if tdecim > 1:
        tmp = total_time_size
        if total_time_size % tdecim == 0:  # necessary because arange doesn't keep the upper limit
            tmp += 1
        time_vector = np.arange(0, tmp, tdecim) * dT + from_time

    # read only one block out of two
    if skip:
        first_record = 0
        last_record = 4 * bts4
        time_indices = range(first_record, last_record, tdecim)

    # read all blocks
    else:
        # starting index (or time sample) in the first block of interest
        first_record = bts4
        # ending index (or last time sample) in the last block of interest
        last_record = 3 * bts4
        time_indices = range(first_record,last_record, tdecim)

    return first_block, last_block,  step_block, time_vector, time_indices


#
# ====================================    GET_SPACE_BOUND()  =======================================
#
def __get_space_bounds__(hd, from_position, to_position, ddecim, section_space):
    '''
    :param hd =  class DasFileHeader instance
    :param from_position:
    :param to_position:
    :param ddecim:
    :param section_space:

    :return: distance_vector selected, array of distance indices in the input distance_vector
    '''

    ZI_start = hd.dist_info['start']
    ZI_end = hd.dist_info['end']
    dX = hd.dist_info['step']

    distance_fiber_in = np.arange(ZI_start, ZI_end + dX, dX)

    if not section_space:
        first_point = 0  # index for 1st point
        last_point = hd.dist_info['npts'] - 1  # index for last point
        distance_fiber_out=distance_fiber_in
        # total_space_size = hd.block_space_size
        # from_position = ZI_start
        # to_position = ZI_end

    else:
        # check range given as argument
        if from_position < distance_fiber_in[0]:
            from_position = distance_fiber_in[0]
        if to_position > distance_fiber_in[-1]:
            to_position = distance_fiber_in[-1]
            print('Warning: drange max is too high, set to ', str(to_position))

        # set indices
        idx = np.where(
            abs(distance_fiber_in - from_position) == min(abs(distance_fiber_in - from_position)))
        first_point = int(idx[0])
        idx = np.where(abs(distance_fiber_in - to_position) == min(abs(distance_fiber_in - to_position)))
        last_point = int(idx[0])
        # total_space_size = (last_point - first_point) + 1


    # distance_length = Total_space_size
    # distance_start = First_point
    # distance_end = distance_start + distance_length
    distance_fiber_out = distance_fiber_in[first_point:last_point + 1:ddecim]  # add 1 because last index is excluded
    distance_indices   = range(first_point, last_point + 1, ddecim)  # add 1 because last index is excluded

    return distance_fiber_out, distance_indices

def __get_time_indices_in_block__(blk, first_blk, last_blk, blk_time_size, first_record, last_record, tdecim):
    '''
    Compute the indices of time samples to be read from current block (or chunck) of data
    We need to distinguish the special cases of the first and last blocks where we read
    starting at first_record and ending at last record
    Things are made a bit tricky when we read with a decimation factor
    :param blk:             current block
    :param first_blk:       first block
    :param last_blk:        last block
    :param blk_time_size:   block size
    :param first_record:    first record to be read in first block
    :param last_record:     last record to be read in last block
    :param tdecim:          decimation factor

    :return: time_indices 
    '''

    bts4 = int(blk_time_size / 4)  # see note in get_time_bounds() to get meaning of bts?
    bts2 = int(blk_time_size / 2)

    # remainder= number of sample to skip at start of a block
    #            when decimating
    # Example: tdecim=3
    # last sample in block N-1 is the last sample of the block => remainder = 2
    # last sample in block N-1 is the one before last sample of the block => remainder = 1
    # last sample in block N-1 is 2 samples before the last => remainder = 0
    #
    if tdecim == 1:
        remainder = 0
    else:
        remainder = int((blk - first_blk) * bts2 - first_record ) % tdecim

    if blk == first_blk:
        time_start = int(first_record)
        time_end = bts2 + bts4

    elif blk == last_blk:
        time_start = bts4
        time_end = bts2 + int(last_record)
    else:
        time_start = bts4
        time_end = bts2 + bts4

    time_start += remainder
    time_indices = range(time_start, time_end, tdecim) #time_end has not to be reached here

    return time_indices

#
# ====================================    INFO()  =======================================
#
def info(filename, verbose=1):
    '''
    Display the informations of a Febus DAS hdf5 file
    :param filename: HDF5 file written by Febus A1 Das interrogator
    :param verbose:  = if 1 (default) only acquisition information is printed; if >1, internal HDF5 file structure is printed
    '''

    Zone_to_Plot = 1

    #
    f = h5py.File(filename, 'r')

    #
    # display groups and datasets
    #
    if verbose > 1:
        print('HDF5 file structure: Groups and Dataset in file ' + filename + ':')
        f.visit(__myprint__)

    # Equivalent de la commande h5disp de matlab
    # ici faire un iterateur sur les cles et dire si c'est un Group ou un Dataset
    # et pour chaque Group, donner les attributs
    #
    if verbose > 1:
        print('\n\nHDF5 file structure, Details :\n------------------------------')
        __display__(f)
        print('')

    #
    # Organisation:    / DAS / SOURCE / ZONE
    #

    # Initialisation
    # Récupération des données d'enetete, dans les attributs du groupe Source1
    #
    print('Acquisition information\n-----------------------')
    DAS = list(f.keys())  # only 1 elemts in this dict.
    DASName = DAS[0]  # 1st key gives name of groupe
    DASNode = f['/' + DASName]
    SOURCES = list(DASNode.keys())
    SOURCEName = SOURCES[0]
    SOURCENode = f['/' + DASName + '/' + SOURCEName]
    SrcAttr = SOURCENode.attrs  # it's a dict.
    ZONEName = 'Zone' + str(Zone_to_Plot)
    ZONENode = f['/' + DASName + '/' + SOURCEName + '/' + ZONEName]
    ZneAttr = ZONENode.attrs

    DAS_name = DASName
    Frequency = SrcAttr['PulseRateFreq'] / 1000.  # search by key
    Freq_res = SrcAttr['FreqRes'] / 1000.  # search by key
    Acquisition_length = SrcAttr['FiberLength']
    Sampling_Res = SrcAttr['SamplingRes']

    Origin = ZneAttr['Origin']
    Spacing = ZneAttr['Spacing']
    Extent = ZneAttr['Extent']
    ZI_start = Origin[0] + Extent[0] * Spacing[0]
    ZI_end = Origin[0] + Extent[1] * Spacing[0]

    SRATEName = '/' + DASName + '/' + SOURCEName + '/' + ZONEName + '/' + 'StrainRate'
    SRATENode = f[SRATEName]
    nb_Block = SRATENode.shape[0]

    #
    print('Device :' + DAS_name)
    print('Fiber Length : ' + str(Acquisition_length) + ' m')
    print('Zone of interest : ' + str(Zone_to_Plot) + ' from ' + str(ZI_start) + ' m to ' + str(ZI_end) + ' m')
    print('Total time of acquisition : ' + str((nb_Block) * Freq_res) + ' s ')
    print('Spatial resolution : ' + str(Spacing[0]) + ' m')
    print('Temporal resolution : ' + str(Spacing[1]) + ' ms')
    print('PRF original : ' + str(Frequency) + ' Hz')
    print('Sampling resolution original : ' + str(Sampling_Res) + ' cm')

    f.close()


def load(filename):
    '''
    Load a binary instance of the DasSection class from disk
    :param filename: (with or without .obj suffix)
    :return: an instance of DasSection class
    '''

    if ".obj" in filename:
        return pickle.load(open(filename, 'rb'))
    else:
        return pickle.load(open(filename + '.obj', 'rb'))

def reduction(filein, fileout, trange=None, drange=None, tdecim=1, ddecim=1, hpcorner=None, kchunk=10, skip=True, verbose=0, filterOMP=False, Zone=1, Source=1):
    '''
    Read a DAS section from an hdf5 Febus-Optics format, extract the requested part or all of it,
    perform a transposition of the section, an optional highpass filter in the space domain
    !!! Decimation are performed without any lowpass filtering !!!
    After reduction, the individual trace can be accessed directly through direct hdf5 requests


    The reducted file is organized as follow
    /
     dataset distance = distance vector
     dataset time =     time vector
     dataset section = section gather of the traces that reference all the traces
     group Traces
        dataset 0 : first trace (single)
        dataset 1 : second trace
        .....

    trace <k> can thus be read either reading directly /Traces/k
    of after reading section as section(k)

    With matlab this looks like:
    ---------------------------
    time     = h5read(file,'/time');      % read the time vector
    distance = h5read(file,'/distance');  % read the distance vector
    tracek   = h5read(file,'/Traces/k');  % read a single trace

    section  = h5read(file,'/section');   % read all traces into a cell array
    tracek   = section{k}

    =======================================================

    :param filein, fileout: hdf5 (.h5) file to read/write
    :param trange:  time range in sec (start,end), (default = None, everything is read)
    :param drange:  distance range in meter (not km) (start,end), (default = None, everything is read)
    :param ddecim:  distance decimation factor (default=1)
    :param tdecim:  time decimation factor (default=1)
    :param hpcorner: Corner frequency for High pass spatial filtering (ex. 600m, default=None)
    :param kchunk:  the ouput HDF5 chunk size is set to input_time_block_size * kchunk (default=10)
    :param skip:    skip redundant block when reading data (default=true)
    :param verbose: be verbose (default=0, minimum message level)
    :param filterOMP: use the multithread sosfilter package instead of scipy.signal.sosfiltfilt
                      (see sosfilter.py to compile fortran source with f2py)
    :param Zone:    not used yet, (default=1)
    :param Source:  not used yet, (default=1)

    '''

    from scipy import signal
    if filterOMP:
        import sosfilter

    #
    #  --------------------------- Time Range option -----------
    #

    if trange:
        section_time = True
        from_time = trange[0]  # in sec
        to_time = trange[1]
    else:
        section_time = False
        from_time = None
        to_time = None

    #
    #  --------------------------- distance Range option -----------
    #

    if drange:
        section_space = True
        from_position = drange[0]
        to_position = drange[1]
    else:
        section_space = False
        from_position = None
        to_position = None

    #
    #  --------------------------- open file for reading
    #
    f = h5py.File(filein, 'r')
    #
    #  --------------------------- open file for writing
    #
    fout = h5py.File(fileout,'w')


    #
    #  --------------------------- read header ----------------------
    #
    hd = DasFileHeader.read(f)

    #
    #  --------------------------- block time size
    #  we read the full block and skip redundant blocks
    #  or read half of it and read all blocks
    if skip:
        block_time_size = int(hd.block_time_size)
    else:
        block_time_size = int(hd.block_time_size / 2)
    #
    # check whether time decimation factor is ok or not, it must divide
    # the chunk size
    #
    if (block_time_size) % tdecim != 0:
        print('Error: time decimation factor must be a divider of the chunk size =',block_time_size)
        print('try :')
        for i in range(2,100):
            if (block_time_size) % i == 0:
                print('  tdecim=',i)
        exit(0)

    #
    #   --------------------------- compute time bounds and time vector -----------------
    #

    # indices for chunk (block) of data
    # indices for first and last time record in the first and last block
    # Vector of times in range[from_time, to_time] with tdecim
    first_block, last_block, step_block, \
    time_out, time_ix = __get_time_bounds__(hd, from_time, to_time, tdecim, section_time, skip)

    #
    #     ---------------------------  compute distance bounds and indices ----------------
    #
    dist_out, dist_ix = __get_space_bounds__(hd, from_position, to_position, ddecim, section_space)

    #
    # ---------------------------   print summary  -------------------
    #
    print(' ')
    print('> Data extraction from [', time_out[0], ' - ', time_out[-1], '] sec and from [', dist_out[0], ' - ',
          dist_out[-1], '] m')
    print('> sampling rate are :',hd.time_info['step'],' sec and ',hd.dist_info['step'],' m')

    #
    # --------------------------- size of data to be written ---------------------------
    #
    total_time_size = len(time_out)
    total_space_size = len(dist_out)

    #
    # ---------------------------- compute optional filter coeffcient
    #
    if hpcorner:
        k_nyquist = np.pi/hd.dist_info['step']
        k_corner = 2*np.pi/hpcorner
        sos = signal.butter(3, k_corner / k_nyquist, 'highpass', output='sos')
        scipy.io.savemat('sos.mat', mdict={'sos': sos})

    #
    #   --------------------------- write header dataset structure on output file
    #   Create one dataset per distance in order to write in a transposed way
    #
    # Each block gives decimed_block_time_size time samples
    # we will write decimed_block_time_size * tdecim * kchunk time sample per chunk
    # A chunk is filled when cblocks are read
    decim_blk_time_size = int(block_time_size/tdecim)
    chunk_size = decim_blk_time_size * tdecim * kchunk
    print('Chunck size is ',chunk_size,' bytes')
    ncblock = kchunk * tdecim

    # create groupe '/Traces'
    grp = fout.create_group('Traces')

    # create dataset that contains reference (!! version >=2.10.0)
    if hasattr(h5py,'ref_dtype'):
        section_ref = fout.create_dataset('section', (total_space_size,), dtype=h5py.ref_dtype)
    else:
        ref_dtype = h5py.special_dtype(ref=h5py.Reference)
        section_ref = fout.create_dataset('section', (total_space_size,), dtype=ref_dtype)

    # create dataset that contains time vector
    time_dset=fout.create_dataset('time',data=time_out,dtype='f8')
    time_dset.attrs.create('dt',hd.time_info['step'])
    # create dataset that contains distance vector
    dist_dset=fout.create_dataset('distance',data=dist_out,dtype='f8')
    dist_dset.attrs.create('dx',hd.dist_info['step'])

    #create datasets that contains traces
    section_list = []
    for i in range(0,total_space_size):
        # define a dataset per spatial location
        dset = grp.create_dataset(str(i), (total_time_size,), chunks=(chunk_size,), dtype='f4')

        # store it in a list
        section_list.append(dset)

        # set dataset attribute
        if hasattr(h5py, 'ref_dtype'):
            dset.attrs.create('H5PATH','/Traces/'+str(i))
        # make link between reference and dataset
        section_ref[i]=dset.ref


    #
    # --------------------------- loop reading blocks ----------------------
    #
    buff_in  = np.empty((decim_blk_time_size, total_space_size), np.float32, 'C')
    buff_out = np.empty((total_space_size, chunk_size), np.float32, 'C')
    #block_list=list(range(first_block, last_block))
    ooffset=0
    i = 0
    for block in range(first_block, last_block, step_block): #in range(0,len(block_list)):

        if verbose >= 1:
            print('    ' + str(block - first_block + 1) + '/' + str(last_block - first_block) + ' blocks', end='\r')

        # copy current block data into buffer
        buff_in[:, :] = hd.srate_node[block, time_ix[0]:time_ix[-1]+1:tdecim,
                                                        dist_ix[0]:dist_ix[-1]+1:ddecim]

        # highpass space filter if requested
        if hpcorner:
            # replace NaN by 0
            buff_in[np.where(np.isnan(buff_in))] = 0.
            if filterOMP:
                # openMP sosfilter version
                buff_in[:, :] = sosfilter.sosfiltfilt(sos, buff_in[:, :])
            else:
                # scipy sequential filter
                buff_in[:, :] = signal.sosfiltfilt(sos, buff_in[:, :], axis=1)


        # transpose data
        buffer_trans = np.transpose(buff_in)

        # Fill output buffer; when it's full, write datasets
        jchunk = i % ncblock
        buff_out[:, jchunk * decim_blk_time_size:(jchunk+1)*decim_blk_time_size] = buffer_trans
        if jchunk == ncblock - 1:
            # write the current time chunk in all spatial dataset
            for j in range(0,total_space_size):
                dset = section_list[j]
                dset[ooffset:ooffset + chunk_size] = buff_out[j,:]
            ooffset += chunk_size
        i += 1

    # end of file reached, write partially filled buffer
    if jchunk != ncblock - 1:
        for j in range(0,total_space_size):
            dset = section_list[j]
            dset[ooffset:ooffset + (jchunk+1)*decim_blk_time_size] = buff_out[j, 0:(jchunk+1)*decim_blk_time_size]
    f.close()
    fout.close()

