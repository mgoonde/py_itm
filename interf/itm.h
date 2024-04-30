typedef struct t_itm_ptr *t_itm_ptr;

t_itm_ptr* itm_create();
void itm_free( t_itm_ptr* itm );
int itm_add_template( t_itm_ptr *, int, int*, double*, int, char*, int );
int itm_match( t_itm_ptr*, int, int*, double*, double, double*, int* );
int itm_set_data( t_itm_ptr*, void*, int, int, int*, void* );
int itm_compute( t_itm_ptr*, double );
int itm_get_dtype( char* );
int itm_get_drank( char* );
int* itm_get_dsize( t_itm_ptr*, char* );
int itm_get_result( t_itm_ptr*, char*, void* );
void itm_print( t_itm_ptr*, int );
void itm_print_canon( t_itm_ptr*, int );
void itm_check_fast( t_itm_ptr* );
char* itm_pg_int2char( int );
void itm_compare_templates( t_itm_ptr*, int, int );
void itm_compare_site_template( t_itm_ptr*, int, int );
void itm_compare_site_canon( t_itm_ptr*, int, int );
void itm_compare_sites( t_itm_ptr*, int, int );



int itm_compute2( t_itm_ptr*, double );
