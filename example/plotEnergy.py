#
#  run test_server.py from another window before running this script
#  O. Coutant, ISTerre, UGA, 2020 from python code by Febus
#  plotting from Test_REQ_DAS_SRPlot.py script by Gaetan Calbris, Febus Optics
#


import a1das
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from scipy.fftpack import fft



def EnergyBand(a,fmin,fmax,nfft):
   
    fs=1/a.data_header['dt']
   
    
    freq=np.linspace(0,fs,nfft)

    a.data = np.nan_to_num(a.data)
    magnitude=abs(fft(a.data,axis=0))
    energy=np.sum(magnitude[np.where((freq>fmin)&(freq<fmax)),:],axis=1)
    return energy
    

# parameters
port = 6668
port = 16667
address = '127.0.0.1'
address = '192.168.1.10'
buffer_Size=50
fmin,fmax=15,50


# open stream
f = a1das.open(a1das.tcp_address(address, port),format='socket')
Distance = f.dist()
Time = f.time()
Size_Time = len(Time)
Size_Dist = len(Distance)

# read 1 block and prepare plot
a = f.read()
nspace = a['nspace']

StrainRate_buff = np.zeros((int(buffer_Size),nspace))

#compute the energy 


E=EnergyBand(a,fmin,fmax,nfft=Size_Time)
StrainRate_buff [0,:] = E

#graphic parameters
E_max = np.max(E)
#print(E)
#print(E.shape)
#E_max=1e7

figure, axes=plt.subplots(figsize=(16,5))
ToPlot = axes.imshow(StrainRate_buff,vmin=0,vmax=E_max ,
                     extent=[min(Distance),max(Distance),0,buffer_Size], aspect= 'auto',
                     origin='lower', cmap='jet')
axes.set_title('Energy')
axes.set_xlabel('Distance [meters]')
axes.set_ylabel('nb_blocks')
figure.colorbar(ToPlot)
plt.show(block=False)

# loop reading further block and plot
Running_Time = 0
Time_Max = 3600*5
nb=0
while Running_Time<Time_Max:
    a = f.read()
    E= EnergyBand(a,fmin,fmax,nfft=Size_Time)
    StrainRate_buff = np.vstack((StrainRate_buff[1:, :],E))
    
    #print(np.max(E))
    ToPlot.set_data(StrainRate_buff)
    figure.canvas.draw_idle()
    figure.canvas.flush_events()
