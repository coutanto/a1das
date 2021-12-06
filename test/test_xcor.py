#%load_ext autoreload
#%autoreload 2
import a1das
from a1das import xcor  #only for debug purpose
from a1das import _a1xcorPy  #only for debug purpose
import numpy as np
import matplotlib.pyplot as plt

# 3 traces au total
ntrace = 4
ntime=1024

# traces = bruit gaussien
traces=np.random.randn(ntime,4)
traces[:,1]=np.roll(traces[:,0],-50)
# On calcule les correlation entre 1 et 2 , 1 et 3
couple=np.array([[1,1],
                 [1,2],
                 [4,1],
                 [1,3]],dtype='int32')


# on enregistre les couples de xcorr à calculer
xcor.register_couple(couple, ntrace, base=1)

print('nbre de correlations ',xcor.get_nxcorr())

## on choisit les options de processing
# help(register_par) donne la liste et les arguments
xcor.register_par('onebit')

# Pour le debug, on regarde la matrice interne de calcul
xcor.infos()
xcor.list_processing()


# On récupère la liste des traces utilisees. cela permet
# de ne sélmectionner qu'une partie des traces pour le calcul des correlations
print('liste des traces utilisées: ',xcor.trace_index_list())


# on calcul les correlations
lag=100
xcor, lags, ijx, ier = xcor.compute_xcorr(traces, lag, verbose=1)

print('ier=',ier)

#
# dessin
#
Time=np.arange(0,ntime)
plt.figure()
plt.plot(Time,traces[:,0],Time,traces[:,1])
plt.figure()
plt.plot(lags,xcor[1,:],lags,xcor[0,:])
plt.show()
