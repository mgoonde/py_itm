cimport numpy as np

cdef extern from "itm.h":
    ctypedef struct t_itm_ptr:
        pass

    t_itm_ptr* itm_create()
    void itm_free( t_itm_ptr* )
    int itm_add_template( t_itm_ptr*, int, int*, double*, int, char*, int )
    int itm_match( t_itm_ptr*, int, int*, double*, double, double*, int* )
    void itm_print( t_itm_ptr* )

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

    def match(self, nat, np.ndarray[int,ndim=1] typ, np.ndarray[double,ndim=2] coords, thr):

        cdef double matched_dh
        cdef int matched_tmplt
        cdef int ierr
        matched_tmplt = itm_match( self.handle, nat, &typ[0], &coords[0,0], thr, &matched_dh, &ierr )
        return matched_tmplt, matched_dh


    def print_all_templates( self ):
        itm_print( self.handle )


    
