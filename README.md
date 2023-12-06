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

   - `cython`:
     ```
     pip3 install cython
     ```
 - The IRA/SOFI library (provided here as git submodule).

## Installation

Clone the repository (SSH/https, you choose):

```
git clone git@github.com:mgoonde/py_itm.git your/path/name
cd your/path/name
```

Initialize and update submodules (only needed the first time after clone):

```
git submodule update --init
```

Compile ITM:

```
make
```

Add the python interface directory to `PYTHONPATH` env variable. Add in your `.bashrc` (edit proper path):

```
export PYTHONPATH=/your/path/to/py_itm/interf:$PYTHONPATH
```

## Run (outdated)

Run the code as:

```
mpiexec -n N python3 py_itm.py -fname FILE -rcut RCUT -shift SHIFT -st SYM_THR -combo DO_COMBO
```

where ALL of the following must be set:
 - `-fname FILE`    : file that ovito can read, i.e. lammps data, xyz, etc.
 - `-rcut RCUT`     : radial cutoff value for local environments (`float` type value).
 - `-shift SHIFT`   : keyword with possible values `gc` or `atm`, it has effect of the found PG, sometimes.
 - `-st SYM_THR`    : threshold for accepted symmetry operations, in units of atomic distance (`float` type value)
 - `-combo DO_COMBO`: logical `T`/`F` whether to perform combinations of the found symmetry operations for each structure.

Example call:

```
mpiexec -n 8 python3 py_itm.py -fname ../run1/dump2.lammpstrj -rcut 2.8 -shift gc -combo T -st 0.3
```


