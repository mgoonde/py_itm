module itm

  use iso_c_binding
  use m_linked_list
  use m_template
  implicit none

  type :: t_itm_ptr
     integer :: ntemplate
     type( linked_list ) :: template_list
   contains
     procedure :: get_template
  end type t_itm_ptr


contains

  function itm_create()result(cptr)bind(C,name="itm_create")
    implicit none
    type( t_itm_ptr ), pointer :: fptr
    type( c_ptr ) :: cptr

    allocate( t_itm_ptr :: fptr )

    !! initialize linked list
    fptr% template_list = linked_list()
    fptr% ntemplate = 0

    cptr = c_loc(fptr)
  end function itm_create

  subroutine itm_free( cptr )bind(C,name="itm_free")
    implicit none
    type( c_ptr ), value :: cptr
    type( t_itm_ptr ), pointer :: fptr

    call c_f_pointer( cptr, fptr )

    call fptr% template_list% destroy()
    deallocate( fptr )
  end subroutine itm_free

  function itm_add_template( cptr, nat_in, typ_in, coords_in, normalize, cmode, cignore_chem ) &
       result(cerr)bind(C, name="itm_add_template")
    implicit none
    type( c_ptr ), value :: cptr
    integer( c_int ), value, intent(in) :: nat_in
    type( c_ptr ), value :: typ_in
    type( c_ptr ), value :: coords_in
    logical( c_bool ), value :: normalize
    type( c_ptr ), value :: cmode
    logical( c_bool ), value :: cignore_chem
    integer( c_int ) :: cerr

    type( t_itm_ptr ), pointer :: fptr
    type( t_template ), pointer :: tmplt
    class(*), pointer :: p

    integer :: nat, i
    integer( c_int ), pointer :: ctyp(:)
    real( c_double ), pointer :: ccoords(:,:)
    integer, dimension(nat_in) :: typ
    real, dimension(3,nat_in) :: coords
    real :: scale
    logical :: rescale, ignore_chem
    character(:), allocatable :: mode

    cerr = 0_c_int
    call c_f_pointer( cptr, fptr )

    !! cast input data
    nat = int( nat_in )
    call c_f_pointer( typ_in, ctyp, shape=[nat_in] )
    call c_f_pointer( coords_in, ccoords, shape=[3,nat_in] )
    allocate( mode, source = c2f_string(cmode) )

    typ = int( ctyp )
    coords = real( ccoords )

    !! check if identical template already exists
    ! do i = 1, fptr% ntemplate
    ! end do


    !! create new template
    nullify( tmplt )
    rescale = logical( normalize )
    ignore_chem = logical( cignore_chem )
    ! write(*,*) "add received ignore_chem:",ignore_chem

    tmplt => t_template( nat, typ, coords, rescale, mode, ignore_chem )

    !! add ptr to list
    fptr% ntemplate = fptr% ntemplate + 1
    nullify(p)
    p => tmplt
    call fptr% template_list% add_pointer( fptr% ntemplate, p )

    deallocate( mode )
  end function itm_add_template


  function itm_get_max_rcut( cptr, cerr )result( rcut )bind(C, name="itm_get_max_rcut")
    implicit none
    type( c_ptr ), value :: cptr
    integer( c_int ), intent(out) :: cerr
    real( c_double ) :: rcut

    integer :: i
    type( t_itm_ptr ), pointer :: fptr
    type( t_template ), pointer :: tmplt
    real :: frcut

    cerr = 0_c_int
    call c_f_pointer( cptr, fptr )

    frcut = -99.9
    !! go through list of templates, get rcut
    do i = 1, fptr% ntemplate
       nullify( tmplt )
       tmplt => fptr% get_template( i )
       !!
       if( .not. associated(tmplt) ) then
          write(*,*) "error getting template pointer",i
          cerr = -1_c_int
          return
       end if

       frcut = max( frcut, tmplt% rcut )
    end do

    !! go through list again, now take the nat and see how far we are

    rcut = real( frcut, c_double )

  end function itm_get_max_rcut


  function itm_match( cptr, nat_in, typ_in, coords_in, thr, matched_dh, cerr )result(matched_tmplt)&
       bind(C,name="itm_match")
    use m_itm_core !, only: natoms_dist_less_than, struc_get_dh, struc_get_scale, struc_rescale
    use m_template
    implicit none
    type( c_ptr ), value :: cptr
    integer( c_int ), value, intent(in) :: nat_in
    type( c_ptr ), value :: typ_in
    type( c_ptr ), value :: coords_in
    real( c_double ), value, intent(in) :: thr
    real( c_double ), intent(out) :: matched_dh
    integer( c_int ), intent(out) :: cerr
    integer( c_int ) :: matched_tmplt

    type( t_itm_ptr ), pointer :: fptr
    type( t_template ), pointer :: tmplt

    integer :: nat, i, candidate_match
    integer( c_int ), pointer :: ctyp(:)
    real( c_double ), pointer :: ccoords(:,:)
    integer, dimension(nat_in) :: typ
    real, dimension(3,nat_in) :: coords

    real :: scale, candidate_match_dh, dh, kmax_factor, fthr
    real, allocatable :: dh_m(:)
    integer :: n, j, hash

    integer, allocatable :: t1(:), t2(:)
    real, allocatable :: c1(:,:), c2(:,:)

    cerr = 0_c_int
    call c_f_pointer( cptr, fptr )

    nat = int( nat_in )
    call c_f_pointer( typ_in, ctyp, shape=[nat] )
    call c_f_pointer( coords_in, ccoords, shape=[3,nat] )

    !! fprec
    fthr = real( thr )
    typ = int( ctyp )
    coords = real( ccoords )

#ifdef DEBUG
    write(*,*) "itm_match received:"
    write(*,*) nat
    write(*,*) "thr",fthr
    do i = 1, nat
       write(*,*) typ(i), coords(:,i)
    end do
#endif


    !! init
    !! array for keeping final dh of matching each struc
    allocate( dh_m(1:fptr% ntemplate), source = 9999.9)
    matched_tmplt = 0
    matched_dh = 999.9

    kmax_factor = 1.3

    !! go through list of templates, match each
    do i = 1, fptr% ntemplate
       nullify( tmplt )
       tmplt => fptr% get_template( i )

       ! write(*,*) "tmplt",i, tmplt% mode
       ! write(*,*) "tmplt% nat", tmplt% nat

       if( .not. associated(tmplt) ) then
          write(*,*) "error getting template pointer",i
          cerr = -1_c_int
          return
       end if

#ifdef DEBUG
       write(*,*) "tmplt:",i
#endif

       !! if n atoms different, cycle this template
       ! if( nat /= tmplt% nat ) cycle


       !! if n equal, first try simple graph isomorphism:
       !! if hash equal, try canon order + svd
       !! if hash not equal, try cshda+svd
       !! if not good, try some common rotations (4 max)

       !! select the subsection of coords to compare, according to template mode
       select case( tmplt% mode )
       case( TEMPLATE_MODE_RCUT )
          !! select n atoms with dist =< rcut
          n = natoms_dist_less_than( nat, coords, tmplt% rcut )

          ! write(*,*) "mode rcut, tmplt% rcut:", tmplt% rcut
          ! write(*,*) "n should be", n

          !!
          !! if this n is diferent from tmplt, cycle
          !! NOTE: eventually this should be modified, because for non-spheric strucs
          !! it will happen practically always. In that case allow ira to match noneq
          if( n /= tmplt% nat ) then
#ifdef DEBUG
             write(*,*) "n different from tmpt% nat",n, tmplt%nat
#endif
             cycle
          end if

       case( TEMPLATE_MODE_NN )
          !! select first N atoms
          n = tmplt% nat

          ! write(*,*) "mode nn, n should be:", n
          !! are there this many atoms in current struc?
          if( nat < n ) then
#ifdef DEBUG
             write(*,*) "nat < n", nat, n
#endif

             cycle
          end if

       case default
          write(*,*) "INVALID Template_mode code:", tmplt% mode
          cerr = -1_c_int
          return

       end select



       !! allocate and set the atomic types
       allocate( t1(1:tmplt% nat), source=tmplt% typ )
       allocate( t2(1:n), source=typ(1:n) )
       if( tmplt% ignore_chem ) then
          !! overwrite all atomic types with equal value
          t1(:) = 1; t2(:) = 1
       end if

       !! allocate and resise coords if needed
       allocate( c1(1:3,1:tmplt% nat), source=tmplt% coords(:,1:tmplt% nat) )
       allocate( c2(1:3,1:n), source=coords(:,1:n) )
       if( tmplt% rescale ) then
          !! get scale of current conf
          scale = struc_get_scale( n, c2 )
          !! rescale the template to match the scale of conf
          call struc_rescale( tmplt% nat, c1, scale )

#ifdef DEBUG
          write(*,*) "tmplt is rescaled"
          write(*,*) tmplt% nat
          write(*,*)
          do j = 1, tmplt% nat
             write(*,*) t1(j), c1(:,j)
          end do
#endif


       end if




       !! Call ira. This should be last resort
       dh = struc_get_dh( &
              tmplt% nat, t1, c1, &
              n, t2, c2, &
              kmax_factor, thr )

#ifdef DEBUG
       write(*,*) "matched dh:", dh, i
#endif

       deallocate( t1, t2 )
       deallocate( c1, c2 )

       dh_m(i) = dh
       !! we matched, no need to go further
       if( dh_m(i) .lt. fthr ) exit

    end do
    ! write(*,*) "minloc dh_m", minloc( dh_m )
    ! write(*,*) "minval dh_m", minval( dh_m )
    candidate_match = minloc( dh_m, 1 )
    candidate_match_dh = minval( dh_m, 1)


    !! the match is below desired thr, set result
    if( candidate_match_dh < thr ) then
       matched_tmplt = candidate_match
       matched_dh = candidate_match_dh
    end if



    deallocate( dh_m )
  end function itm_match


  subroutine itm_print(cptr)bind(C,name="itm_print")
    !! print all templates in the list
    implicit none
    type( c_ptr ), value :: cptr
    type( t_itm_ptr ), pointer :: fptr

    type( t_template ), pointer :: t
    integer :: i

    call c_f_pointer( cptr, fptr )

    do i = 1, fptr% ntemplate
       t => fptr% get_template( i )
       write(*,*) "got node",i
       call t% print()
    end do

  end subroutine itm_print


  !!-----------------------------------------------------------
  !! local functions

  function get_template( fptr, idx )result( template_pointer )
    !! get template :: get pointer to node in template_list
    implicit none
    class( t_itm_ptr ), intent(in) :: fptr
    integer, intent(in) :: idx
    type( t_template ), pointer :: template_pointer

    class(*), pointer :: p

    nullify( template_pointer, p )

    call fptr% template_list% get( idx, p)
    if( .not.associated(p) ) then
       write(*,*) "get_template::p not associated, node does not exist?"
       return
    end if

    select type( p )
    class is( t_template )
       template_pointer => p
    class default
       write(*,*) "get_template::type p is wrong?"
       return
    end select

  end function get_template

  FUNCTION c2f_string(ptr) RESULT(f_string)
    use iso_c_binding

    INTERFACE
       !! standard c function
       FUNCTION c_strlen(str) BIND(C, name='strlen')
         IMPORT :: c_ptr, c_size_t
         IMPLICIT NONE
         TYPE(c_ptr), INTENT(IN), VALUE :: str
         INTEGER(c_size_t) :: c_strlen
       END FUNCTION c_strlen
    END INTERFACE

    TYPE(c_ptr), INTENT(IN) :: ptr
    CHARACTER(LEN=:), ALLOCATABLE :: f_string
    CHARACTER(LEN=1, KIND=c_char), DIMENSION(:), POINTER :: c_string
    INTEGER :: n, i

    IF (.NOT. C_ASSOCIATED(ptr)) THEN
       f_string = ' '
    ELSE
       n = INT(c_strlen(ptr), KIND=KIND(n))
       CALL C_F_POINTER(ptr, c_string, [n+1])
       allocate( CHARACTER(LEN=n)::f_string)
       do i = 1, n
          f_string(i:i) = c_string(i)
       end do
    END IF
  END FUNCTION c2f_string

end module itm
