import sysconfig
import numpy as np
import os

#python include directory
dic = sysconfig.get_paths()
print('PYTHON_INC = ',dic['include'])

#binary modules suffix
suffix=sysconfig.get_config_var('SOABI')
print('SUFFIX = ',suffix+'.so')

inc= np.get_include()
print('NUMPY_INC = ',inc)

#os.system('f2py --help-link f2py_info | grep sources | sed -e "s/sources/F2PY_SRC/" | sed -e "s/\[//" | sed -e "s/\]//"')
#os.system('f2py --help-link f2py_info | grep include_dirs | sed -e "s/include_dirs/F2PY_INC/" | sed -e "s/\[//" | sed -e "s/\]//"')
