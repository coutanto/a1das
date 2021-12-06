#
# Example of a client reading real-time data from DAS or from a file
# and plotting them as a movie
#
#
#  Plot traces as curves (vector version)
#  O. Coutant, ISTerre, UGA, 2020 
#  read <nblock> at a time and plot without overlap
#  plot <ntrace> max
#  Variable area option is usually too slow and loose data


import a1das
import matplotlib.pyplot as plt

# parameters
# for the DAS
port = 16667
address = '192.168.1.10'
# for a file, set according to the parameter chosen when running a1das.vserver.ZMQVirtualDasServer()
# usually, port = 6667 and address = 127.0.0.1

va=False
nblock=4
clip=50
ntrace=50

# open stream
f = a1das.open(a1das.tcp_address(address, port),format='socket')

# read 2 block and make 1st plot
a = f.read(block=nblock)
fig, lines = a.plot(max=ntrace,clip=clip, variable_area=va)

Running_Time = 0
Time_Max = 3600*5
while Running_Time<Time_Max:
    a = f.read(block=nblock)
    a.plot(fig=fig, max=ntrace, redraw=lines, clip=clip, variable_area=va)
