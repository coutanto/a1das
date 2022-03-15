__doc__="Various tools to manipulate das data"


def save_location_as_h5(filename, position, dist, trace, crs, pos_label, tgte=None):
    """
    ## Description
    Save (or append) fiber location information to a h5 file by adding/replacing the group "/location"
    and adding/replacing attributes.

        content of the group /location:
        datasets:
        /location/position (x,y,z)
        /location/trace    trace_index
        /location/dist     curvilinear abscissa
        attributes:
        crs: coordinate reference system
        pos_label: (x,y,z) labels
        tgte: optional, tangent vector in 3D space

    ## Input
        filename: (string) file to store location information. If data_filename is not defined, it is assumed to
                  be the same as filename, and location information is appended to it.
        position: (ndarray, float) [3 x ntraces] (x,y,z) or (lat,lon,z) or (.,.,.) for each trace
        dist:     (ndarray, float) [ntraces] curvilinear abscissa
        trace:   (ndarray, int)  [ntraces] index of traces
        crs: (string) coordinate reference  system using EPSG ids (see https://epsg.io, ex "EPSG:3857" for
             WGS 84 / Pseudo-Mercator)
        labels: (tuple) lables for the three coordinates (ex: ('x','y','z') or ('lat','lon','elev')
    """

    import h5py
    from .core import open
    from ._a1das_exception import WrongValueError, DataTypeError
    from numpy import ndarray

    #
    # check args
    #
    if not isinstance(position, ndarray):
        raise DataTypeError('positions must be a ndarray')
    if position.shape[0] != 3:
        raise DataTypeError('positions must be a ndarray of shape (3,ntrace)')
    if not isinstance(trace, ndarray):
        raise DataTypeError('itrace must be a ndarray')
    if not isinstance(crs, str):
        raise DataTypeError('crs (coordinate reference system) must be a string')
    if not isinstance(pos_label, tuple) or len(pos_label) != 3:
        raise DataTypeError('labels must be a tuple of 3 strings')
    for i,label in enumerate(pos_label):
        if not isinstance(label,str):
            raise DataTypeError('labels must be a tuple of strings')

    npos = position.shape[1]
    ntrace = trace.shape[0]
    ndist = dist.shape[0]
    if npos != ntrace or npos != ndist or ndist != ntrace:
        raise DataTypeError('ndist, ntrace and npos must have the same size')
    #
    # add/replace position dataset
    # dataset is created as resizable with unlimited length
    #
    try:
        f = h5py.File(filename, 'r+')
        g = f['/location']
        loc_attr = g.attrs
        loc_dset = g['position']
        loc_dset.resize((3,npos))
        dist_dset = g['distance']
        dist_dset.resize((npos,))
        trace_dset = g['trace']
        trace_dset.resize((npos,))
        if tgte is not None:
            try:
                tgte_dset = g['tgte']
                tgte_dset.resize((3,npos))
                tgte_dset = tgte
            except:
                tgte_dset = g.create_dataset("tgte", data=tgte, maxshape=(3, None))

        loc_dset[::] = position
        dist_dset[:] = dist
        trace_dset[:] = trace

    except:
        f = h5py.File(filename, 'w')
        g = f.create_group("location")
        loc_dset = g.create_dataset("position", data=position, maxshape=(3,None))
        dist_dset = g.create_dataset("distance", data=dist, maxshape=(None,))
        trace_dset = g.create_dataset("trace", data=trace, maxshape=(None,), dtype='int32')
        if tgte is not None:
            tgte_dset = g.create_dataset("tgte", data=tgte, maxshape=(3,None))
        loc_attr = g.attrs

    loc_attr['crs'] = crs
    loc_attr['pos_label'] = pos_label

    f.close()

def read_location_from_h5(filename):
    """
    ## Description
    Read fiber location information from a h5 file

    ## Return
        positions: ndarray of size 3 x ntrace (x,y,z for each trace)
        dist: ndarray of size ntrace          (curvilinear distance for each trace)
        trace: ndarray of size ntrace         (trace index)
        tgte: (optional) ndarray of size 3 x ntrace (tangent vector to the fiber)
        crs= (str) reference coordinate system
        pos_label= (tuple of str) labels for the 3 positions ('x','y','z') or ('lat','lon','elev') or ...
    """
    import h5py
    from ._a1das_exception import FileFormatError

    f = h5py.File(filename,'r')

    try:
        g = f['/location']
        crs = g.attrs['crs']
        pos_label = g.attrs['pos_label']
        position = g['position'][::]
        dist = g['distance'][:]
        trace = g['trace'][:]
        try:
            tgte = g['tgte'][::]
        except:
            tgte = None
    except:
        raise FileFormatError('file '+filename+'does not contain a location group')

    f.close()
    if tgte is None:
        return position, dist, trace, crs, pos_label
    else:
        return position, dist, trace, tgte, crs, pos_label

def save_location_as_wkt(fileout, positions, crs, distance=None):
    """
    ## Definition
    Save location information in a .csv file using WellKnownText convention to be read by QGIS.
    Positions are supposed to be passed as (x,y,z) or (lat,lon,z) according to reference coordinate system

    ## Input
        fileout= (str) filename without suffix
        positions= (ndarray) of size 3 x ntrace
        rcs= (str) reference coordinate system in the form EPGS:#### or ####, only the number #### is used
        distance= distance along fiber (optional)
    """
    from ._a1das_exception import WrongValueError

    # check dimension
    ntrace1 = positions.shape[1]
    try:
        ntrace2 = distance.shape[0]
        if ntrace1 != ntrace2:
            raise WrongValueError('position and distance size are not consistent')
    except:
        pass

    ix = crs.find(':') + 1
    srid = crs[ix:]
    with open(fileout+'.csv','w') as f:
        f.write('SRID='+str(srid)+'\n')
        f.write('WKT,fid\n')
        if distance is not None:
            for i in range(0,ntrace1):
                print('\"POINT ZM (%.2f %.2f %.2f %.2f)\",\"%d\"' % (positions[0,i],positions[1,i],positions[2,i],distance[i],i),file=f)
        else:
            for i in range(0,ntrace1):
                print('\"POINT Z (%.2f %.2f %.2f)\",\"%d\"' % (positions[0,i],positions[1,i],positions[2,i],i), file=f)


def save_location_as_gpkg(fileout, position, dist, trace, crs, tgte = None):
    """
    ## Description
    Save fiber location information using the geo-package format (.gpkg) which can be directly imported into QGIS.
    This requires the 'geopandas' package
    ## Input
    position: (ndarray, float) [3 x nspace] x,y,z where (x,y) depends on the reference coordinate system
    dist:  (ndarray, float) [nspace] curvilinear abscissa
    trace: (ndarray, int) [nspace] trace index
    crs: (str) coordinate reference system in the standard form 'epgs:####'
    tgte: optional,(ndarray, float) [3 x nspace]; tangent vector to the fiber
    """
    import geopandas as gpd
    from shapely.geometry import point
    from numpy import ndarray
    from ._a1das_exception import DataTypeError

    if ((not isinstance(position, ndarray)) or (not isinstance(dist, ndarray)) or
        not isinstance(trace, ndarray)):
        raise DataTypeError('expect numpy.ndarray as input arguments')
    if not isinstance(crs, str):
        raise DataTypeError('expect Coordinate reference system to be a string of the form "epgs:#####')

    npos = position.shape[1]
    ntrace = trace.shape[0]
    ndist = dist.shape[0]
    if npos != ntrace or npos != ndist or ndist != ntrace:
        raise DataTypeError('ndist, ntrace and npos must have the same size')

    filename=[ None for d in dist]
    geometry = gpd.points_from_xy(position[0,:], position[1,:],position[2,:])
    if tgte is None:
        gdf=gpd.GeoDataFrame(data={'curvilinear_abscissa':dist, 'trace':trace, 'filename':filename,
                                   'geometry':geometry}, crs = crs)
    else:
        gdf=gpd.GeoDataFrame(data={'curvilinear_abscissa':dist, 'trace':trace,
                                   'tgtx':tgte[0,:],'tgty':tgte[1,:],'tgtz':tgte[2,:],
                                   'filename':filename,
                                   'geometry':geometry}, crs = crs)
    gdf.to_file(fileout+'.gpkg', driver='GPKG')

def read_location_from_gpkg(file):
    """
    ##Description
    Read fiber location information from a .gpkg (geopackage) file

    ##Return
    position: a numpy ndarray of size [3 x npoints] or [2 x npoints] if elevation is not known
    crs: coordinate reference system
    """
    import geopandas as gpd
    from shapely.geometry import Point
    from numpy import array,asarray
    from ._a1das_exception import DataTypeError

    gdf = gpd.read_file(file)
    # retrieve only the first geometry object
    path_geom = gdf.geometry[0]
    # check it's a collection of points
    if not isinstance(path_geom, Point):
        raise DataTypeError('wrong geometry property in .gpkg file, expect Point')
    else:
        # path_geom is a Shapely Point object (see Shapely python module)
        # and coords is the field that allows to access the points that defines
        # the line
        #path_xyz = asarray(path_geom).T
        #or
        path_xyz = asarray([g.coords for g in  gdf.geometry]).squeeze().T
        return path_xyz, gdf.crs.srs


def dt_morlet(sig1, sig2, wlen, freq):
    """
    ## Description
    Compute time shift between sig1 and sig2 using a moving projection on morlet wavelet transform

    ## Input
    ----------
    sig1 : ndarray, first signal input
    sig2 : ndarray, second signal input
    wlen : int, window length.
    freq : float or 3 elements float list, frequency bandwidth of morlet ex: f; or [f1,f2] or [fmin, fmax, df]

    ## Returns
    dt: ndarray(nfreq, ntime), float, time shift vector (=phase / omega)
    ws1, ws2: ndarray(nfreq, ntime), complex, projection on Morlet wavelet
    """
    from scipy.signal import morlet, correlate
    from numpy import arange, asarray, pi, angle, ndarray, unwrap, dot, sqrt
    
    #
    # set frequency range
    #
    if isinstance (freq, float):
        freq = asarray([freq])
    if (len(freq)==2):
        freq = asarray(freq)
    if len(freq)==3:
        freq = arange(freq[0],freq[1]+freq[2],freq[2])

    #
    # initialize  output array
    #
    dt = ndarray((len(freq),len(sig1)))
    ws1 = ndarray((len(freq),len(sig1)),dtype=complex)
    ws2 = ndarray((len(freq),len(sig2)),dtype=complex)
    for i,f in enumerate(freq):
        omega = 2*pi*f
        wm = morlet(wlen, omega)
        norm =  sqrt(sum(abs(wm)**2))
        wm = wm/norm
        ws1[i,:] = correlate(sig1, wm, mode='same')
        ws2[i,:] = correlate(sig2, wm, mode='same')
        dt[i,:] = unwrap(angle(ws2[i,:])-angle(ws1[i,:]))/omega

        
    return freq, dt, ws1, ws2
        

