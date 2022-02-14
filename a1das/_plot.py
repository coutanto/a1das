#
# Module to plot sections
# O. Coutant, ISTerre, UGA, 2020
#
#
# Content:
#
#         plot() = plot a section as vectors
#        rplot() = plot section as raster
#

import matplotlib.pyplot as plt

__doc__='Functions calld by core.plot and core.rplot'
#
# ====================================    RPLOT()  =======================================
#
def rplot(a1, fig=None, cmap='RdBu', vrange=None, splot=(1, 1, 1), title='', drange=None, trange=None):
    '''
    plot a DAS section as a raster image
    ------------------------------------
    a1: an A1Section class instance
    fig = a fig handle necessary to plot several plot on the same figure (default=None, a new figure is created)
    cmap = colormap (default='RdBu')
    vrange = (vmin,vmax) range for colormap (default=[] set to max value)
    splot = subplot position (default=(1,1,1))
    '''
    import numpy as np

    if not fig:
        fig = plt.figure(figsize=(18, 9))

    transpose = False
    if 'data_axis' in a1.data_header.keys(): #old format
        if a1.data_header['data_axis'] == 'space_x_time':
            transpose = True
    if 'axis1' in a1.data_header.keys():
        if a1.data_header['axis1'] == 'space':
            transpose = True

    dist = a1.data_header['dist']
    time = a1.data_header['time']
    if drange is None:
        d1=0
        if transpose:
            d2 = a1.data.shape[0]
        else:
            d2=a1.data.shape[1]
        dmin=dist[0]
        dmax=dist[-1]
    else:
        dmin=drange[0]
        dmax=drange[1]
        d1=np.nanargmin(np.abs(dist-dmin))
        d2=np.nanargmin(np.abs(dist-dmax))

    if trange is None:
        t1=0
        if transpose:
            t2 = a1.data.shape[1]
        else:
            t2 = a1.data.shape[0]
        tmin=time[0]
        tmax=time[-1]
    else:
        tmin=trange[0]
        tmax=trange[1]
        t1=np.nanargmin(np.abs(time-tmin))
        t2=np.nanargmin(np.abs(time-tmax))

    #check for nan
    a1.data = np.nan_to_num(a1.data)

    if not vrange:
        if transpose:
            vmax = np.max(np.absolute(a1.data[d1:d2, t1:t2]), axis=None)
        else:
            vmax = np.max(np.absolute(a1.data[t1:t2, d1:d2]), axis=None)
        vmin = -vmax
    else:
        vmin = vrange[0]
        vmax = vrange[1]

    ax = plt.subplot(splot[0], splot[1], splot[2])

    if transpose:
        extent = [tmin, tmax, dmax, dmin]
        pos = ax.imshow(a1.data[d1:d2, t1:t2], cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto', extent=extent)
        ax.set_ylabel('distance meters')
        ax.set_xlabel('time sec')
    else:
        extent = [dmin, dmax, tmax, tmin]
        pos = ax.imshow(a1.data[t1:t2,d1:d2], cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto', extent=extent)
        ax.set_xlabel('distance meters')
        ax.set_ylabel('time sec')

    ax.set_title('section DAS ' + title)
    fig.colorbar(pos)

    return fig, ax
#
# ====================================    PLOT()  =======================================
#
def plot(a1, fig=None, clip=100, splot=(1, 1, 1), title='', max=100, amax=None, by_trace=False, variable_area=False,
         drange=None, trange=None, redraw=None):
    """
    plot a DAS section as curves
    ----------------------------
    a1: an A1Section class instance
    fig   = a fig handle necessary to plot several plot on the same figure (default=None, a new figure is created)
    clip  = % of max amplitude that is clipped (deault=100) (clip and amax are exclusice)
    amax  = maximum amplitude that is clipped (default=None)
    splot = subplot position (default=(1,1,1))
    max=  = max number of trace per plot along distance dimension (default=100)
    by_trace = (bool) normalize by plot (False) or by trace (True) (default=False)
    variable_area = variable area plot, positive wiggles are filled (default=False)
    drange= distance range
    trange= time range
    redraw = {lines} redraw lines of a precedent figure  or None for new plot, requires fig argument
    """

    import numpy as np
    import matplotlib.dates as mdates
    from datetime import datetime, timezone

    if fig is None:
        fig = plt.figure(figsize=(18, 9))

    if redraw is not None:
        lines = redraw
        redraw = True
        ax = fig.axes[0]
    else:
        ax = plt.subplot(splot[0], splot[1], splot[2])
        redraw = False
        lines={}

    dist = a1['dist']
    time = a1['time'] + a1['otime']
    otime, otimes = a1.otime()

    transpose = False
    if 'data_axis' in a1.data_header.keys(): #old format
        if a1.data_header['data_axis'] == 'space_x_time':
            transpose = True
    if 'axis1' in a1.data_header.keys():
        if a1.data_header['axis1'] == 'space':
            transpose = True


    if drange is None:
        d1=0
        if transpose:
            d2 = a1.data.shape[0]
        else:
            d2=a1.data.shape[1]
        dmin=dist[0]
        dmax=dist[-1]
    else:
        dmin=drange[0]
        dmax=drange[1]
        d1=np.nanargmin(np.abs(dist-dmin))
        d2=np.nanargmin(np.abs(dist-dmax))
        if d2==d1:
            d2=d1+1

    if trange is None:
        t1=0
        if transpose:
            t2 = a1.data.shape[1]
        else:
            t2=a1.data.shape[0]

    else:
        tmin=trange[0]
        tmax=trange[1]
        t1=np.nanargmin(np.abs(time-tmin - a1['otime']))
        t2=np.nanargmin(np.abs(time-tmax - a1['otime']))

    from_distance = dmin
    to_distance = dmax

    # Do not plot more than MAX traces
    # adjust decimation rate on distance
    #ndist = len(a1.dist)  # number of distance
    ndist = d2-d1+1
    if ndist < max:
        ddecim = 1
    else:
        ddecim = int(np.fix(ndist / max))

    #check for nan
    a1.data = np.nan_to_num(a1.data)

    #idrange = list(range(0, ndist, ddecim))
    idrange = list(range(d1, d2, ddecim))
    ntraces = len(idrange)
    if transpose:
        max_value = np.max(np.absolute(a1.data[d1:d2, t1:t2]), axis=None)
    else:
        max_value = np.max(np.absolute(a1.data[t1:t2, d1:d2]), axis=None)
    clip = clip / 100
    if not amax:
        clipped_value = max_value * clip
    else:
        clipped_value = amax


    width = 2 * clipped_value * ntraces
    if width == 0.:
        width = 1.
    offset = np.arange(0, ntraces + 1) * 2 * clipped_value
    gain = (to_distance - from_distance) / width
    if gain == 0:
        gain = 1

    # compute the width, trace offset, and gain
    dates =[datetime.fromtimestamp(ts, tz=timezone.utc) for ts in time[t1:t2]]
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S.%f'))
    for i in range(0, len(idrange)):
        if transpose:
            v = a1.data[idrange[i], t1:t2].copy()
        else:
            v = a1.data[t1:t2, idrange[i]].copy()
        if by_trace:
            max_val = np.max(v)
            v = v /max_val*max_value
            #if not amax:
            #    clipped_value = max_val * clip
            #else:
            #    clipped_value = amax
        v[v > clipped_value] = clipped_value
        v[v < -clipped_value] = -clipped_value
        v = v + offset[i]
        v *= gain
        v += from_distance
        if not redraw:
            lines[i], = ax.plot(dates,v,'k')
        else:
            lines[i].set_ydata(v)
            lines[i].set_xdata(dates)
            ax.set_xlim(dates[0],dates[-1])

        if variable_area:
            trace_offset = gain*offset[i]+from_distance
            ax.fill_between(dates, trace_offset, v, where=(v > trace_offset), color='k')
        del v

    if redraw:
        fig.canvas.draw()
        fig.canvas.flush_events()
    else:
        plt.show(block=False)

    if not redraw:
        ax.set_xlabel('time sec')
        ax.set_ylabel('distance m')
        ax.set_title('section DAS, start_time= ' + otimes + title)


    return fig, lines

