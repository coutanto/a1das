#import a1das_util
import a1das
import matplotlib.pyplot as plt



#fichier reduit
filei='SR_2021-01-14_11-55-54_UTC.h5'
fileo='SRRNT_2021-01-14_11-55-54_UTC.h5'
a1das.reduction.reduction_notranspose(filei,fileo)
fileo='SRRT_2021-01-14_11-55-54_UTC.h5'
a1das.reduction.reduction_transpose(filei,fileo)

f1=a1das.open('SRRNT_2021-01-14_11-55-54_UTC.h5','reducted')
f2=a1das.open('SRRT_2021-01-14_11-55-54_UTC.h5','reducted')

a1=f1.read()
a2=f2.read()

a1.plot()
a1.rplot()
a2.plot()
a2.rplot()

plt.show()
