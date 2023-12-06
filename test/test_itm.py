import numpy as np
import py_itm
from ovito.io import import_file
from ovito.data import CutoffNeighborFinder
from mpi4py import MPI
import argparse
from func import *


# init MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# parser
p = argparse.ArgumentParser()

p.add_argument( "-fname", "--filename", help="input filename", required=True, type=str )
p.add_argument( "-rcut", help="radial cutoff", required=True, type=float )
p.add_argument( "-dthr", help="distance threshold", required=False, type=float, default=0.3 )
args = p.parse_args()

filename = args.filename
rcut = args.rcut
dthr = args.dthr



# init itm, different instance on each rank
itm=py_itm.itm()

# read template struc from file
nat, typ, coords, options = read_xyz( "cubd2_2n" )
# add template to all mpi ranks
ierr = itm.add_template( nat, typ, coords, options )

nat, typ, coords, options = read_xyz( "int_2n" )
ierr = itm.add_template( nat, typ, coords, options )

# nat, typ, coords, options = read_xyz( "int2_4.4" )
# ierr = itm.add_template( nat, typ, coords, options )

nat, typ, coords, options = read_xyz( "vac1_2n" )
ierr = itm.add_template( nat, typ, coords, options )


# if rank ==0:
#     itm.print_all_templates()

comm.barrier()

# import sys
# sys.exit()

if rank == 0:
    # itm.print_all_templates()

    ## Read file using ovito function
    f=import_file( filename )

    ## create DataCollection object for ovito (maybe does nothing, not sure)
    dc=f.compute()
    nat = dc.particles.count
    # nat = 5
    typ = dc.particles.particle_types
    coords = dc.particles.position
    lat = dc.cell

    #
    ## set the finder function
    #

    # obtain max rcut from templates
    ## NOTE: cannot do when any template is rescaled, default to -rcut from input
    # rcut1 = itm.get_max_rcut()
    # add 10%
    rcut *= 1.1
    # print( "max rcut:", rcut)

    # define neighbor finder by this rcut, sort by distances
    neigh = CutoffNeighborFinder( rcut, dc ).find_all( sort_by='distance' )

    neiglist=neigh[0]
    veclist=neigh[1]

    ## cast data into regular int/float
    nat = int(nat)
    typ = np.cast[int](typ)
    neiglist = np.cast[int](neiglist)
    veclist = np.cast[float](veclist)
    coords = np.cast[float](coords)

    print(nat, flush=True )
    print( 'Lattice="',lat[0][0], lat[0][1], lat[0][2], lat[1][0], lat[1][1], lat[1][2], lat[2][0], lat[2][1], lat[2][2], '" properties=species:I:1:pos:R:3:id:I:1:structure_type:I:1:dHausdorff:R:1', flush=True)

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
# coords = comm.bcast( coords, root=0)

# number of structures to compute per cpu rank
per_rank = int(nat/(size))
# maximal index reached with the split to cores:
# print( "M",(size)*per_rank)
# there can be modulo index uncomputed, if the number of atoms is indivisible by mpi size
rest = nat - size*per_rank
# print("rest:",rest)

matched_arr = np.zeros([nat], dtype=int)
matched_dhs = np.zeros([nat], dtype=float)

## Are we off-by-one, OBIWAN? Kennobe, we are correct.
ob1 = 0


# loop over slice of indices per_rank
for i in range( ob1+rank*per_rank, ob1+(rank+1)*per_rank):
# for i in range( 1, 5 ) :

    ## read from neigh list, also gives coords relative to idx=i
    idx, coords_loc = extract_elements( neiglist, veclist, i )
    nat_loc = len(idx)
    typ_loc = np.array(typ[idx],dtype=np.int32)
    if nat_loc == 2:
        # this shoudl return error
        print("site:",i, "has too few neighbors:",nat_loc, flush=True)
        msg = "Error, local env too small, see message above."
        raise ValueError( msg )

    matched_template, matched_dh = itm.match( nat_loc, typ_loc, coords_loc, dthr )

    # v = coords[i]
    # print( "%i %14.8f %14.8f %14.8f %i %s" %(typ[i], v[0], v[1], v[2], i, pg ) )
    # print( "%i %14.8f %14.8f %14.8f %i %i %14.10f" %(typ[i], v[0], v[1], v[2], i, matched_template, matched_dh ) )
    # print( "%i %14.8f %14.8f %14.8f %i %i %14.10f" %(typ[i], v[0], v[1], v[2], i, matched_template, matched_dh ), flush=True )

    matched_arr[i] = matched_template
    matched_dhs[i] = matched_dh
    # print( "%i %i %i" %(i, typ[i], matched_template ) )

    # print("====",flush=True)
    # print( "m",nat_loc, flush=True)
    # print("m",matched_template,flush=True)
    # for j,v in enumerate( coords_loc):
    #     print("m", typ_loc[j], v[0], v[1], v[2], flush=True )
    # print("====",flush=True)

# if the last few atoms have been missed, do them now
if rank == 0:
    if rest > 0:
        for i in range( nat-1, nat-1-rest, -1 ):
            ## read from neigh list, also gives coords relative to idx=i
            idx, coords_loc = extract_elements( neiglist, veclist, i )
            nat_loc = len(idx)
            typ_loc = np.array(typ[idx],dtype=np.int32)
            if nat_loc == 2:
                # this shoudl return error
                print("site:",i, "has too few neighbors:",nat_loc, flush=True)
                msg = "Error, local env too small, see message above."
                raise ValueError( msg )
            matched_template, matched_dh = itm.match( nat_loc, typ_loc, coords_loc, dthr )
            # v = coords[i]
            # print( "%i %14.8f %14.8f %14.8f %i %s" %(typ[i], v[0], v[1], v[2], i, pg ) )
            # print( "%i %14.8f %14.8f %14.8f %i %i %14.10f" %(typ[i], v[0], v[1], v[2], i, matched_template, matched_dh ) )
            matched_arr[i] = matched_template
            matched_dhs[i] = matched_dh


tot_i = np.zeros( nat, dtype=int )
comm.Allreduce( [matched_arr, MPI.INT], [tot_i, MPI.INT], op=MPI.SUM)
tot_d = np.zeros( nat, dtype=float)
comm.Allreduce( [matched_dhs, MPI.DOUBLE], [tot_d, MPI.DOUBLE], op=MPI.SUM)

comm.barrier()
if rank ==0:
    for i, v in enumerate(coords):
        print( "%i %14.8f %14.8f %14.8f %i %s %11.6f" %(typ[i], v[0], v[1], v[2], i, tot_i[i], tot_d[i] ) )
