This directory contains the `cython` extenstion module which acts as interface to call the
functions from `libitm.so`.

Files needed for the compilation:
  - `setup.py` : definitions where to find the lib etc.;
  - `itm.h` : C-header for the functions inside `libitm.so`.
  - `py_itm.pyx` : the extension module which declares and calls functions from `libitm.so`;

Manual compilation by:
```
python3 setup.py build_ext --inplace
```


Import the module into python as:
```
 ## python3
 import py_itm
 me = py_itm.itm()
```
To make the module visible from any directory, edit the `PYTHONPATH`:
```
export PYTHONPATH=/path/to/py_itm/interf:$PYTHONPATH
```

NOTE:
Why `cython` and not `ctypes`? Supposedly `cython` is faster, but I don't see a big difference.
