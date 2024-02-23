
import py_itm
itm=py_itm.itm()

from ovito.data import CutoffNeighborFinder
from ovito.io import import_file
import sys
f=import_file( sys.argv[1] )
dc=f.compute()
n=CutoffNeighborFinder( 3.5, dc ).find_all( sort_by="distance" )

itm.set_data( dc, n )
itm.compute(dthr=0.5)
itm.check_fast()

t = itm.get_result( "site_template" )
dh = itm.get_result( "site_dh" )
pg = itm.get_result( "site_pg" )

lat=dc.cell.matrix



with open("tt.xyz", "w") as f:
   print( dc.particles.count, file=f )
   print('Lattice="%.6f %.1f %.1f %.1f %.6f %.1f %.1f %.1f %.6f " properties=species:I:1:pos:R:3:id:I:1:structure_type:I:1:dHausdorff:R:1:molecule_type:I:1' %(lat[0][0], lat[0][1], lat[0][2], lat[1][0], lat[1][1], lat[1][2], lat[2][0], lat[2][1], lat[2][2]), file=f )
   for i, v in enumerate( dc.particles.positions ):
       print( "%i %.6f %.6f %.6f %i %i %.6f %s"%(dc.particles.particle_types.array[i], \
               v[0], v[1], v[2], i+1, t[i], dh[i], pg[i]), file=f )


