import a1das
import h5py
import numpy as np

filei='Raw_2021-02-03_09-56-57_UTC.h5'
fileo='strain_2021-02-03_09-56-57_UTC.h5'
a1das.reduction.raw2strain(filei, fileo, 2., 0.005, order_space=2,use_compression=True,kchunk=1,verbose=2, transpose=False)
fileo='strain__transposed_2021-02-03_09-56-57_UTC.h5'
a1das.reduction.raw2strain(filei, fileo, 2., 0.005, order_space=2,use_compression=True,kchunk=1,verbose=2, transpose=True)
