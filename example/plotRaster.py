#
# Example of a client reading real-time data from DAS or from a file
# and plotting them as a movie
#
#
#  run test_server.py from another window before running this script
#  O. Coutant, ISTerre, UGA, 2020 from python code by Febus
#  plotting from Test_REQ_DAS_SRPlot.py script by Gaetan Calbris, Febus Optics
#


import a1das
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# parameters
# for the DAS
port = 16667
address = '192.168.1.10'
# for a file, set according to the parameter chosen when running a1das.vserver.ZMQVirtualDasServer()
# usually, port = 6667 and address = 127.0.0.1
# parameters

strain_max = 5000
buffer_Size=10

# open stream
f = a1das.open(a1das.tcp_address(address, port),format='socket')
Distance = f.dist()
Time = f.time()
Size_Time = len(Time)
Size_Dist = len(Distance)

# read 1 block and prepare plot
a = f.read(block=[])

StrainRate_buff = np.zeros((int(buffer_Size*Size_Time),int(Size_Dist)))
StrainRate_buff [:Size_Time,:] = a.data

#graphic parameters
figure, axes=plt.subplots(1)
ToPlot = axes.imshow(StrainRate_buff,vmin=-strain_max,vmax=strain_max,
                     extent=[min(Distance),max(Distance),min(Time),max(Time*buffer_Size)], aspect= 'auto',
                     origin='lower', cmap=cm.seismic)
axes.set_title('StrainRate (nm/m/s)')
axes.set_xlabel('Distance (m)')
axes.set_ylabel('Time (s)')
figure.colorbar(ToPlot)
plt.show(block=False)

# loop reading further block and plot
Running_Time = 0
Time_Max = 3600*5
nb=0
while Running_Time<Time_Max:
    a = f.read()
    StrainRate_buff = np.vstack((StrainRate_buff[Size_Time:, :], a.data))
    ToPlot.set_data(StrainRate_buff)
    figure.canvas.draw_idle()
    figure.canvas.flush_events()
#    plt.pause(0.1)
