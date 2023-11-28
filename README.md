# py_itm

Characterizing local atomic environments using IRA/SOFI.
Playing with MPI in python3.

## Requirements

 - The following python3 modules:
   - `ovito`:

     ```
       pip3 install ovito
     ```

   - `mpi4py`:
     ```
       pip3 install mpi4py
     ```

 - The IRA/SOFI library (provided here as git submodule).

## Installation

Clone the repository (SSH/https, you choose):

```
    git clone git@github.com:mgoonde/py_itm.git your/path/name;
    cd your/path/name
```

Initialize and update submodules (only needed the first time after clone):

```
    git submodule update --init
```

Compile IRA/SOFI:

```
    cd IRA/src ; make all; cd -
```

Add IRA python interface directory to `PYTHONPATH` env variable. Add in your `.bashrc`:

```
    export PYTHONPATH=/your/path/name/IRA/interface:$PYTHONPATH
```
