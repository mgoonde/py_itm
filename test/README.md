# Test of pyitm:

at the beginning, define 3 template structures which are read from files: `cubd_2n`, `int_2n`, and `vac1_2n`. Then match these templates over all the input structure.

### Run
Run test as:
```
mpiexec -n 2 python3 test_itm.py -fname conf.inp -rcut 4.0
```

### Arguments
The argument options are:
 - `-fname`: `str`, required, input filename containing the configuration, any file type that ovito can read.
 - `-rcut`: `float`, required, the radial cutoff distance, which should be such that all structure templates that you define, should be possible to find with the given rcut value. Larger values are generally ok, but will make the code run slower.
 - `-dthr`: `float`, optional (default `dthr = 0.3`), distance threshold at which two structures are considered equal (in the sense of maximal displacement of any one atom, *i.e.* the Hausdorff distance).


### Template files
The template files are written in the `xyz` format, where the second line contains the options for that template.
Possible options for a template are:
 - `-mode`: `str`, optional (default `mode=nn`), the mode used to build the local configuration for this template. Possible values `nn` and `rcut`.
 - `-no_normalize`: `bool`, optional (default `no_normalize=False`), flag to skip the rescaling of templates. This is automatically turned on for any template with `mode=rcut`.
 - `-ignore_chem`: `bool`, optional (default `ignore_chem=False`), flag to ignore the chemical (atomic) types of that template. If `True`, the geometry of the structure will attempt to be matched, disregarding the atomic types.


### Output
The output is an extended-xyz format, readable by ovito. The output properties are `StructureType` and `dHausdorff`.
The `StructureType` property has the value of the index of the template which that local configurations matches to, *e.g.* `StructureType=2` means the atom has been matched to the 2nd template that has been added in the beginning of the `test_itm.py` program. The value of `dHausdorff` is the value of the Hausdorff distance with which this atom has matched to given template.

