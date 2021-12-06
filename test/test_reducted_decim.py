import a1das
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal


#
# create reducted files
#
filei='SR_2021-01-14_11-55-54_UTC.h5'

fileo='SRRT_2021-01-14_11-55-54_UTC.h5'
a1das.reduction.reduction_transpose(filei,fileo)

tdecim=5
fileo='SRRTD_2021-01-14_11-55-54_UTC.h5'
a1das.reduction.reduction_transpose(filei,fileo,tdecim=tdecim)

#
# read and compare
#
f1=a1das.open('SRRT_2021-01-14_11-55-54_UTC.h5','reducted')
f2=a1das.open('SRRTD_2021-01-14_11-55-54_UTC.h5','reducted')

a1=f1.read()
a2=f2.read()

S = a1.obspy_stream(drange=[60,60], station_name='DAS1')
S += a2.obspy_stream(drange=[60,60], station_name='DAS2')

S.plot()

plt.show()

#
# compare spectra for trace at dist=60m
#
id = a1.index(60.)[0]
dt = a1['dt']
n = a1['ntime']
A1 = np.fft.fft(a1.data[id],norm='ortho')#/np.sqrt(n)
freq1 = np.fft.fftfreq(n, d=dt)

dt = a2['dt']
n = a2['ntime']
A2 = np.fft.fft(a2.data[id],norm='ortho')#/np.sqrt(n)
freq2 = np.fft.fftfreq(n, d=dt)
plt.figure()
plt.semilogy(freq1,np.abs(A1),freq2,np.abs(A2))
plt.show()

#
#compare low freq part
#
f_nyquist  = 1./a1['dt'] / 2.
f_corner  = 0.7 * f_nyquist / tdecim
sos = signal.butter(6, f_corner / f_nyquist, 'lowpass', output='sos')
tracef = signal.sosfiltfilt(sos, a1.data[id])

plt.figure()
plt.plot(a1['time'],tracef,a1['time'][0:-1:tdecim],tracef[0:-1:tdecim],a2['time'],a2.data[id])
plt.legend(['undecimated lowpass trace','decimated lowpass trace','reducted decimated trace']) 
plt.show()


