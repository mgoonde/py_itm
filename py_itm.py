from ovito.io import import_file
from ovito.data import CutoffNeighborFinder
import ira_mod
import numpy as np
import argparse
from mpi4py import MPI

# local file
from func import *

# init sofi
sofi = ira_mod.SOFI()

# init MPI
comm=MPI.COMM_WORLD
size=comm.Get_size()
rank=comm.Get_rank()

# parser
p = argparse.ArgumentParser()

# -fname FILENAME -sh SHIFT_TYPE -st SYM_THR -rcut RCUT -combo GET_COMBOS
p.add_argument("-fname","--filename",help="Input filename")
p.add_argument("-shift","--shift_type",help="Type of shift: gc/atm")
p.add_argument("-st","--sym_thr",help="Symmetry threshold for SOFI")
p.add_argument("-rcut","--rcut",help="local env cutoff")
p.add_argument("-combo","--get_combos",help="make combos of found Ops")

# parse arguments
args = p.parse_args()

filename = args.filename
rcut = float( args.rcut )
sym_thr = float( args.sym_thr )
shift = args.shift_type
get_combos = args.get_combos


if rank == 0:
    ## Read file using ovito function
    f=import_file( filename )

    ## create DataCollection object for ovito (maybe does nothing, not sure)
    dc=f.compute()
    nat = dc.particles.count
    # nat=50
    typ = dc.particles.particle_types
    coords = dc.particles.position
    lat = dc.cell

    ## set the finder function
    neigh = CutoffNeighborFinder( rcut, dc ).find_all( sort_by='index' )

    neiglist=neigh[0]
    veclist=neigh[1]

    ## cast data into regular int/float
    nat = int(nat)
    typ = np.cast[int](typ)
    neiglist = np.cast[int](neiglist)
    veclist = np.cast[float](veclist)
    coords = np.cast[float](coords)

    print(nat, flush=True )
    print( 'Lattice="',lat[0][0], lat[0][1], lat[0][2], lat[1][0], lat[1][1], lat[1][2], lat[2][0], lat[2][1], lat[2][2], '" properties=species:I:1:pos:R:3:id:I:1:molecule_type:I:1', flush=True)

else:
    nat=0
    typ=None
    neiglist=None
    veclist=None
    coords = None

comm.barrier()

# bcast the data
nat = comm.bcast(nat)
neiglist = comm.bcast( neiglist, root=0)
veclist = comm.bcast( veclist, root=0)
typ = comm.bcast( typ, root=0)
coords = comm.bcast( coords, root=0)

# define output array
# ppg = np.ndarray([nat], dtype="U4")

# number of structures to compute per cpu rank
per_rank = int(nat/(size))


# if rank ==0:
#     for i in range( ob1+rank*per_rank, ob1+(rank+1)*per_rank):
#         print(i)
# print()
# comm.Barrier()
# if rank ==1:
#     for i in range( ob1+rank*per_rank, ob1+(rank+1)*per_rank):
#         print(i)

## Are we off-by-one, OBIWAN? Kennobe, we are correct.
ob1 = 0

# loop over slice of indices per_rank
for i in range( ob1+rank*per_rank, ob1+(rank+1)*per_rank):
    ## read from neigh list, also gives coords relative to idx=i
    idx, coords_loc = extract_elements( neiglist, veclist, i )
    nat_loc = len(idx)
    typ_loc = typ[idx]
    if nat_loc == 2:
        # this shoudl return error
        print("site:",i, "has too few neighbors:",nat_loc, flush=True)
        msg = "Error, local env too small, see message above."
        raise ValueError( msg )

    # type of shift
    if shift == 'gc':
        # geometric center
        gc = np.mean( coords_loc )
    elif shift =='atm':
        # atom at center (this is already the case -> zero shift)
        gc = np.array([0.0, 0.0, 0.0])

    # shift coords
    coords_loc = coords_loc - gc

    # compute default (this automatically includes combos)
    # sym = sofi.compute( nat_loc, typ_loc, coords_loc, 0.2 )

    nm, m = sofi.get_symm_ops( nat_loc, typ_loc, coords_loc, sym_thr )

    # do combos?
    if get_combos:
       nm, m = sofi.mat_combos(nm, m)

    # find pg name
    pg = sofi.get_pg( nm, m )
    # ppg[i] = pg

    v = coords[i]
    print( "%i %14.8f %14.8f %14.8f %i %s" %(typ[i], v[0], v[1], v[2], i, pg ) )


# if rank == 0:
#     comm.receive( )
# else:
#     comm.send()

# tot = np.ndarray( nat, dtype = "U4" )

# # allreduce with op=SUM is ok
# comm.Allreduce( [ppg, MPI.CHAR], [tot, MPI.CHAR], op=MPI.SUM)

# comm.Barrier()
# if rank ==0:
#     # print( "rank 0 has:")
#     # for i in range(nat):
#     #     print(i, tot[i])

#     print(nat )
#     print( 'Lattice="',lat[0][0], lat[0][1], lat[0][2], lat[1][0], lat[1][1], lat[1][2], lat[2][0], lat[2][1], lat[2][2], '" properties=species:I:1:pos:R:3:id:I:1:molecule_type:I:1')
#     for i, v in enumerate( coords ):
#         print( typ[i], v[0], v[1], v[2], i, tot[i])

# comm.Disconnect()
