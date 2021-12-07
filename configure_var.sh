
#This shell try to determine the good locations and the suffix for applying f2py


python configure_var.py > Makefile.python

f2py --help-link f2py_info | grep sources | sed -e "s/sources/F2PY_SRC/" | sed -e "s/\[//" | sed -e "s/\]//" >> Makefile.python
f2py --help-link f2py_info | grep include_dirs | sed -e "s/include_dirs/F2PY_INC/" | sed -e "s/\[//" | sed -e "s/\]//" >> Makefile.python
