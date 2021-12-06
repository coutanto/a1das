#
# Example of a client reading real-time data from DAS or from a file
# and plotting them as a movie
#
#  Plot traces as curves (vector version)
#  O. Coutant, ISTerre, UGA, 2020 
#  read <nblock> at a time and plot with a 50% overlap
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
va=False   #Variable area plot, slow !!!
nblock=2   # read nblock at a time, 1 block is usually = 1 sec
ntrace=100 # plot at most ntrace
clip=50 #% of max amp
drange=[2000,2500] # distance range

# open stream
f = a1das.open(a1das.tcp_address(address, port),format='socket')

# read 2 block and make 1st plot
a = f.read(block=nblock)
b = f.read(block=nblock)
c = a.concat(b)
fig, lines = c.plot(max=ntrace,drange=drange,clip=clip, variable_area=va)
a = b

Running_Time = 0
Time_Max = 3600*5
while Running_Time<Time_Max:
    b = f.read(block=nblock) # read one more block
    c = a.concat(b)    # concatenate
    a = b
    c.plot(fig=fig, max=ntrace, drange=drange,redraw=lines, clip=clip, variable_area=va)
