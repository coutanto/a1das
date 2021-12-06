#
# simulate a DAS TCP stream 
#

import a1das

# run server for SR data
# file to send block per block
file = '/Volumes/terrain/taconnaz/das/SR_2021-08-26_11-55-53_UTC.h5'
a1das.vserver.ZMQVirtualDasServer(file)

#run server for RAW data
#file = '/media/coutant/DAS-1/data/taconnaz/das/Raw_2021-08-27_09-39-28_UTC.h5'
#gl=4.
#dt=0.005
#a1das.vserver.ZMQRawVirtualDasServer(file, gl, dt, drange=[20,75])

