typedef struct t_itm_ptr *t_itm_ptr;

t_itm_ptr* itm_create();
void itm_free( t_itm_ptr* itm );
int itm_add_template( t_itm_ptr *, int, int*, double*, int, char*, int );
int itm_match( t_itm_ptr*, int, int*, double*, double, double*, int* );
void itm_print( t_itm_ptr* );

