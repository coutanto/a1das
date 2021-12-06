#import a1das_util
import a1das
import matplotlib.pyplot as plt
import numpy as np

filename='SR_2021-01-14_11-55-54_UTC.h5'

# open febus file
f = a1das.open(filename,format='febus')

#test reading febus file
a1=f.read()

#test index function
dlist3 = a1.index(drange=[[dist] for dist in np.arange(50,60.,1.)])

#test setting origin time
a1.set_otime_from_filename()
print(a1)

#test conversion in obspy
S=a1.obspy_stream(drange=[60,60])
print('\n\nnext steps may take a while...large plots coming soon')
S.plot()

#test a1das plot
plt.figure()
a1.plot(trange=[5.,7.])
a1.rplot()
print('\n\nkill figures to keep going')
plt.show()

print('lus des données ',a1.data.shape)

plt.figure()
plt.subplot(3,1,1)
plt.plot(a1.time(),a1.data[:,10])


a2=f.read(trange=[10.2, 21.7], skip=False)


plt.subplot(3,1,2)
plt.plot(a2.time(),a2.data[:,10])

a3=f.read(trange=[10.2, 21.7])


plt.subplot(3,1,3)
plt.plot(a2.time(),a2.data[:,10],a3.time(),a3.data[:,10])
plt.show()

f.close()

#
# test save 
#
a1.save('test_save.h5')
