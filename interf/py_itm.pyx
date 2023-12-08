cimport numpy as np
import numpy as np
from cython.view cimport array as cyarray
from cpython cimport array

cdef extern from "itm.h":
    ctypedef struct t_itm_ptr:
        pass

    t_itm_ptr* itm_create()
    void itm_free( t_itm_ptr* )
    int itm_add_template( t_itm_ptr*, int, int*, double*, int, char*, int )
    int itm_match( t_itm_ptr*, int, int*, double*, double, double*, int* )
    int itm_set_data( t_itm_ptr*, char*, int, int, int*, void* )
    int itm_compute( t_itm_ptr*, double )
    int itm_get_dtype( char* )
    int itm_get_drank( char* )
    int* itm_get_dsize( t_itm_ptr*, char* )
    int itm_get_result( t_itm_ptr*, char*, void* )
    void itm_print( t_itm_ptr* )
    void free(void*)

cdef class itm:
    cdef t_itm_ptr* handle

    def __cinit__(self):
        self.handle = itm_create()

    def close(self):
        if self.handle is not NULL:
            itm_free( self.handle )
            self.handle = NULL

    def __dealloc__(self):
        self.close()

    def add_template( self, nat, np.ndarray[int,ndim=1] typ, np.ndarray[double,ndim=2] coords,\
                      template_options ):
        import argparse
        p = argparse.ArgumentParser()
        p.add_argument( '-no_normalize', default=False, action='store_true')
        p.add_argument( '-mode', type=str, default="nn", choices =["rcut", "nn"] )
        p.add_argument( '-ignore_chem', default=False, action = 'store_true')

        args = p.parse_args( template_options )
        # print(args)
        norm = not args.no_normalize
        mode = args.mode
        ignore_chem = args.ignore_chem

        cdef int ierr
        ierr = itm_add_template( self.handle, nat, &typ[0], &coords[0,0], norm, mode.encode(), ignore_chem )
        return ierr

    def set_data( self, dc, neigh ):
        ## take the things we need from DataCollection,
        ## and the neigh object from ovito

        cdef np.ndarray[int, ndim=1]csize
        cdef int cerr
        cdef int ctyp
        cdef int crank

        ## number of atoms
        ## is in dc.particles.count
        name = "nat"
        cdef int ival = np.intc( dc.particles.count )
        csize = np.array([0], dtype=np.int32)
        ctyp=itm_get_dtype( name.encode() )
        crank=itm_get_drank( name.encode() )
        cerr = itm_set_data( self.handle, name.encode(), ctyp, crank, &csize[0], &ival )

        ## atomic types
        ## is in dc.particles.particle_types
        name = "typ"
        cdef np.ndarray[int,ndim=1]i1d = np.intc(dc.particles.particle_types.array)
        csize = np.array([dc.particles.count], dtype=np.int32)
        ctyp=itm_get_dtype( name.encode() )
        crank=itm_get_drank( name.encode() )
        cerr = itm_set_data( self.handle, name.encode(), ctyp, crank, &csize[0], &i1d[0] )

        ## neigh list
        ## is in neigh[0]
        name = "neighlist"
        cdef np.ndarray[int,ndim=2]i2d = np.intc( neigh[0] )
        csize=np.array([ np.size(neigh[0], 1), np.size(neigh[0], 0)], dtype=np.int32)
        ctyp=itm_get_dtype( name.encode() )
        crank=itm_get_drank( name.encode() )
        cerr = itm_set_data( self.handle, name.encode(), ctyp, crank, &csize[0], &i2d[0,0] )

        ## veclist
        ## is in neigh[1]
        name = "veclist"
        cdef np.ndarray[double,ndim=2]r2d = np.float64( neigh[1] )
        csize=np.array([ np.size(neigh[1], 1), np.size(neigh[1], 0)], dtype=np.int32)
        ctyp=itm_get_dtype( name.encode() )
        crank=itm_get_drank( name.encode() )
        cerr = itm_set_data( self.handle, name.encode(), ctyp, crank, &csize[0], &r2d[0,0] )


    def set_data1( self, name, object val ):
        cdef int cerr=-1

        cdef int ival
        cdef np.ndarray[int,ndim=1]i1d
        cdef np.ndarray[int,ndim=2]i2d
        cdef np.ndarray[double,ndim=2]r2d

        cdef np.ndarray[int, ndim=1] csize
        cdef int crank
        cdef int ctyp

        if name == "nat":
            ival = np.intc(val)
            ctyp = 1
            crank = 0
            csize=np.array([0], dtype=np.int32)
            cerr = itm_set_data( self.handle, name.encode(), ctyp, crank, &csize[0], &ival )
        elif name == "typ":
            i1d = np.intc(val)
            ctyp = 1
            crank = 1
            csize = np.array([np.size(val)], dtype = np.int32)
            if np.ndim( val ) != crank:
                msg = "wrong rank in input! Expected: 1, given:",np.ndim(val)
                raise ValueError( msg )
            cerr = itm_set_data( self.handle, name.encode(), ctyp, crank, &csize[0], &i1d[0] )

        elif name == "neighlist":
            i2d = np.intc(val)
            ctyp = 1
            crank = 2
            csize=np.array([ np.size(val, 1), np.size(val, 0)], dtype=np.int32)
            if np.ndim( val ) != crank:
                msg = "wrong rank in input! Expected: 2, given:",np.ndim(val)
                raise ValueError( msg )
            cerr = itm_set_data( self.handle, name.encode(), ctyp, crank, &csize[0], &i2d[0,0] )
        elif name == "veclist":
            r2d = np.float64( val )
            ctyp = 2
            crank = 2
            csize = np.array([np.size(val,1), np.size(val,0)], dtype=np.int32)
            if np.ndim( val ) != crank:
                msg = "wrong rank in input! Expected: 2, given:",np.ndim(val)
                raise ValueError( msg )
            cerr = itm_set_data( self.handle, name.encode(), ctyp, crank, &csize[0], &r2d[0,0] )


        if cerr != 0:
            msg="error in set_data"
            raise ValueError( msg )
        return


    def compute( self, dthr=None ):
        # default dthr=0.3
        if dthr is None:
            dthr = 0.3
        cdef double cthr=dthr
        cdef int cerr
        cerr = itm_compute( self.handle, cthr )

    def match(self, nat, np.ndarray[int,ndim=1] typ, np.ndarray[double,ndim=2] coords, thr):

        cdef double matched_dh
        cdef int matched_tmplt
        cdef int ierr
        matched_tmplt = itm_match( self.handle, nat, &typ[0], &coords[0,0], thr, &matched_dh, &ierr )
        return matched_tmplt, matched_dh

    def get_dtype( self, name ):
        cdef int dtype
        dtype = itm_get_dtype( name.encode() )
        return dtype

    def get_drank( self, name ):
        cdef int drank
        drank = itm_get_drank( name.encode() )
        return drank

    def get_dsize( self, name ):
        cdef int drank = itm_get_drank( name.encode() )
        cdef int* csize = itm_get_dsize( self.handle, name.encode() )
        psize=np.ndarray([drank],dtype=int)
        for i in range(drank):
            psize[i] = csize[i]
        return psize

    def get_result( self, name ):
        cdef int cerr=0

        # cdef int ctyp =itm_get_dtype( name.encode() )
        # cdef int crank=itm_get_drank( name.encode() )
        # cdef int* csize=itm_get_dsize( self.handle, name.encode() )
        ctyp = self.get_dtype( name )
        crank = self.get_drank( name )
        psize = self.get_dsize( name )

        cdef void* cval
        cerr = itm_get_result( self.handle, name.encode(), &cval )

        cdef int* int_arr
        # cdef np.ndarray[np.int32_t, ndim=1] ipyval = np.empty(csize[0], dtype=np.int32)
        cdef np.ndarray[np.int32_t, ndim=1] ipyval = np.empty(psize, dtype=np.int32)
        if ctyp == 1:
            int_arr = <int*>cval
            # ipyval[:] = <int[:csize[0]]>int_arr
            ipyval[:] = <int[:psize[0]]>int_arr
            free(cval)
            return ipyval

        cdef double* double_arr
        # cdef np.ndarray[np.float64_t, ndim=1] dpyval = np.empty(csize[0], dtype=np.float64)
        cdef np.ndarray[np.float64_t, ndim=1] dpyval = np.empty(psize, dtype=np.float64)
        if ctyp == 2:
            double_arr = <double*>cval
            # dpyval[:] = <double[:csize[0]]>double_arr
            dpyval[:] = <double[:psize[0]]>double_arr
            free(cval)
            return dpyval

    def print_all_templates( self ):
        itm_print( self.handle )


    
