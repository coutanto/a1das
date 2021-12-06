#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 16:36:56 2021

@author: coutant UGA 
"""
#%%
import a1das
import obspy
import obspy.signal as sig
import matplotlib.pyplot as plt
import numpy as np
from a1das import xcor
from time import process_time

#
# comm parameters
#
#For the DAS
port = 16667
address = '192.168.1.10'
#for a file, set according to the server launched by a1das.vserver.ZMQ????
port = 6667
address='127.0.0.1'

ref_unit_is_dist = True
nblock=5 #1 block is 1sec
lag=1.5


# open stream
f = a1das.open(a1das.tcp_address(address, port),format='socket')
# shortcut for data header
dhd = f.data_header
ntrace = dhd['nspace']
ntime = dhd['ntime']
dist = dhd['dist']
print('ntrace = ',ntrace)
print('ntime = ',ntime)

#  definition par indice
site_ref = [500]#, 20, 40, 60, 80]#, 120, 130]
site_target = range(300,700)

#  definition par distance
dist_ref = 700.  # in meter
dist_target = dist #all distance availables
dist_target=np.arange(200,dist[-1]-200,5)

if ref_unit_is_dist:
    tmp = dhd.index([dist_ref, dist_ref])
    site_ref = [tmp[0]]
    site_target = dhd.index(list(dist_target))
    
# set xcorr couple                    
couple =  [list([site_ref[0],j]) for j in site_target]

for iref in site_ref[1:]:
    #site_target = 
    couple += [list([1,j]) for j in site_target]
    
#couple = np.array(couple)
xcor.register_couple(couple,ntrace=ntrace)
# For debug
#xcor.infos()
#xcor.list_processing()



figure, axes=plt.subplots(1)
plt.show(block=False)

Running_Time = 0
Time_Max = 3600*5
init=True
while Running_Time<Time_Max:
    a1 = f.read(block=nblock)
    otime, otimed = a1.otime()
    print('compute xcorr at time ',otimed)
    xcorr, lags, ijx, ier = xcor.compute_xcorr(a1,lag=lag,stack=True,verbose=0)
    xcor.register_couple(couple,ntrace=ntrace)

    # Loop for each xcorr sources
    for i0 in site_ref:
        xcorrs = xcor.sort_xcorr(i0, ijx, xcorr)
        xmax = np.max(np.abs(xcorrs))/10
        if ref_unit_is_dist:
            extent = [lags[0], lags[-1], dist_target[-1], dist_target[0]]
        else:
            extent = [lags[0], lags[-1], trace_cible[-1], trace_cible[0]]
        if init:
            ToPlot = axes.imshow(xcorrs,cmap='seismic',vmin=-xmax,vmax=xmax,aspect='auto',extent=extent)
            if ref_unit_is_dist:
                axes.plot(lags,0.*lags+dist_ref,'r')
                axes.set_ylabel('traces distance (m)')
            else:
                axes.plot(lags,0.*lags+i0,'r')
                axes.set_ylabel('traces index')
            axes.set_title('xcorr '+otimed)
            axes.set_xlabel('lag (s)')
            figure.colorbar(ToPlot)
            plt.show(block=False)
            init = False
        else:
            ToPlot.set_data(xcorrs)
            axes.set_title('xcorr '+otimed)
            figure.canvas.draw_idle()
            figure.canvas.flush_events()
