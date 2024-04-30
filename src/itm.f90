module itm

  use iso_c_binding
  use m_linked_list
  use m_datatypes
  use m_template
  use m_canon
  use m_itm
  implicit none

contains

  function itm_create()result(cptr)bind(C,name="itm_create")
    implicit none
    type( t_itm_ptr ), pointer :: fptr
    type( c_ptr ) :: cptr

    fptr => itm_init()
    cptr = c_loc(fptr)
  end function itm_create

  subroutine itm_free( cptr )bind(C,name="itm_free")
    implicit none
    type( c_ptr ), value :: cptr
    type( t_itm_ptr ), pointer :: fptr

    call c_f_pointer( cptr, fptr )

    call fptr% template_list% destroy()
    call fptr% fast_list% destroy()
    if( allocated(fptr% veclist))deallocate( fptr% veclist )
    if( allocated(fptr% neighlist))deallocate( fptr% neighlist)
    if( allocated( fptr% typ))deallocate( fptr% typ )
    if( allocated( fptr% count_n))deallocate( fptr% count_n )
    if( allocated( fptr% site_template))deallocate( fptr% site_template )
    if( allocated( fptr% site_dh))deallocate( fptr% site_dh )
    if( allocated( fptr% site_pg))deallocate( fptr% site_pg )
    ! if( allocated( fptr% site_strain))deallocate( fptr% site_strain )
    deallocate( fptr )
  end subroutine itm_free

  !! itm_add_template C-wrapper
  function citm_add_template( cptr, nat_in, typ_in, coords_in, normalize, cmode, cignore_chem ) &
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

    integer :: nat
    integer( c_int ), pointer :: ctyp(:)
    real( c_double ), pointer :: ccoords(:,:)
    integer, dimension(nat_in) :: typ
    real, dimension(3,nat_in) :: coords
    logical :: rescale, ignore_chem
    character(:), allocatable :: mode

    call c_f_pointer( cptr, fptr )

    !! cast input data
    nat = int( nat_in )
    call c_f_pointer( typ_in, ctyp, shape=[nat_in] )
    call c_f_pointer( coords_in, ccoords, shape=[3,nat_in] )
    allocate( mode, source = c2f_string(cmode) )

    typ = int( ctyp )
    coords = real( ccoords )
    rescale = logical( normalize )
    ignore_chem = logical( cignore_chem )

    cerr = int( itm_add_template( fptr, nat, typ, coords, rescale, mode, ignore_chem ), c_int )

    deallocate( mode )
  end function citm_add_template


  function itm_set_data( cptr, cname, ctyp, crank, csize, cval)result(cerr)bind(C, name="itm_set_data")
    use m_itm_core, only: count_nn
    implicit none
    type( c_ptr ), value :: cptr
    type( c_ptr ), value :: cname
    integer( c_int ), value :: ctyp
    integer( c_int ), value :: crank
    type( c_ptr ), value :: csize
    type( c_ptr ), value :: cval
    integer( c_int ) :: cerr

    type( t_itm_ptr ), pointer :: fptr
    integer :: dtyp, drank
    character(:), allocatable :: fname
    integer( c_int ), pointer :: isize(:)
    integer( c_int ), pointer :: iptr, i1d(:), i2d(:,:)
    real( c_double ), pointer :: r1d(:), r2d(:,:)

    integer :: i

    cerr = -1_c_int
    call c_f_pointer( cptr, fptr )

    allocate( fname, source=c2f_string(cname) )
    dtyp = int( ctyp )
    drank = int( crank )

    ! write(*,*) "got name:",fname
    ! write(*,*) dtyp, drank



    select case( fname )
    case( "nat" )
       !! check if dtyp is as expected
       if( dtyp /= ITM_DTYPE_INT ) return
       !! check if rank is as expected
       if( drank /= 0 ) return
       !! get data
       call c_f_pointer( cval, iptr )
       fptr% nat = int( iptr )
       !!
       !! allocate arrays for results, if exist destroy
       !!
       if( allocated( fptr% site_template))deallocate( fptr% site_template)
       allocate( fptr% site_template(1:fptr% nat), source=0)
       !!
       if( allocated( fptr% site_dh))deallocate( fptr% site_template)
       allocate( fptr% site_dh(1:fptr% nat), source = 99.9)
       !!
       ! if( allocated( fptr% site_strain))deallocate( fptr% site_strain)
       ! allocate( fptr% site_strain(1:fptr% nat), source=0.0)
       !!
       if( allocated( fptr% site_pg)) deallocate( fptr% site_pg)
       allocate( fptr% site_pg(1:fptr% nat),source=0 )
       !!
       if( allocated( fptr% perm_site2rough))deallocate(fptr% perm_site2rough)
       allocate( fptr% perm_site2rough(1:fptr% nat))

    case( "typ" )
       !! check if dtyp is as expected
       if( dtyp /= ITM_DTYPE_INT ) return
       !! check if rank is as expected
       if( drank /= 1 ) return
       !! get size
       call c_f_pointer( csize, isize, shape=[drank] )
       !! get data
       call c_f_pointer( cval, i1d, shape=isize )
       !! overwrite previous data
       if( allocated(fptr% typ)) deallocate( fptr% typ )
       allocate( fptr% typ, source=int(i1d) )

    case( "neighlist" )
       !! check if dtyp is as expected
       if( dtyp /= ITM_DTYPE_INT ) return
       !! check if rank is as expected
       if( drank /= 2 ) return
       !! get size
       call c_f_pointer( csize, isize, shape=[drank] )
       !! get data
       call c_f_pointer( cval, i2d, shape=isize )
       !! overwrite previous data
       if( allocated( fptr% neighlist))deallocate( fptr% neighlist)
       allocate( fptr% neighlist, source=int(i2d) )
       !! switch to 1-based index
       fptr% neighlist = fptr% neighlist + 1
       ! do i = 1, 40
       !    write(*,*) fptr% neighlist(:,i)
       ! end do
       !! get the count
       allocate( fptr% count_n(1:fptr% nat), source=0 )
       call count_nn( size( fptr% neighlist, 1), size(fptr% neighlist, 2), &
            fptr% neighlist, fptr% nat, fptr% count_n )


    case( "veclist" )
       !! check if dtyp is as expected
       if( dtyp /= ITM_DTYPE_REAL ) return
       !! check if rank is as expected
       if( drank /= 2 ) return
       !! get size
       call c_f_pointer( csize, isize, shape=[drank] )
       !! get data
       call c_f_pointer( cval, r2d, shape=isize )
       !! overwrite previous data
       if( allocated( fptr% veclist))deallocate( fptr% veclist)
       allocate( fptr% veclist, source=real(r2d))
       ! do i = 1, 40
       !    write(*,*) fptr% veclist(:,i)
       ! end do


    case default
       write(*,*) "name not recognized: ", fname
       return
    end select


    cerr = 0_c_int
    deallocate( fname )
  end function itm_set_data


  function itm_compute( cptr, cthr )result( cerr )bind(C, name="itm_compute" )
    use m_itm_core
    use omp_lib
    implicit none
    type( c_ptr ), value :: cptr
    real( c_double ), value :: cthr
    integer( c_int ) :: cerr

    type( t_itm_ptr ), pointer :: fptr
    integer :: i, j, me, t
    integer :: nthread

    integer :: nat_loc, n
    integer, allocatable :: idx_loc(:)
    integer, allocatable :: typ_loc(:)
    real, allocatable :: coords_loc(:,:)

    type( t_template ), pointer :: tmplt

    integer, allocatable :: t1(:), t2(:)
    real, allocatable :: c1(:,:), c2(:,:)
    real :: scale, dh, dh_old, dthr
    integer, allocatable :: site_template(:)
    real, allocatable :: site_dh(:)

    call c_f_pointer( cptr, fptr )

    dthr = real( cthr )
    fptr% dthr = dthr

    !$OMP PARALLEL
    nthread = OMP_GET_NUM_THREADS()
    !$OMP END PARALLEL

    ! write(*,*) "got threads",nthread
    cerr = -1_c_int
    !! check if all data is set
    ! write(*,*) "nat =",fptr% nat
    ! write(*,*) "typ",allocated( fptr% typ )
    ! write(*,*) "veclist",allocated( fptr% veclist )
    ! write(*,*) "neighlist",allocated( fptr% neighlist )

    fptr% site_dh(:) = 99.9
    fptr% site_template(:) = 0
    fptr% site_pg(:) = 0

    !$OMP PARALLEL PRIVATE(me, site_template, site_dh )
    me = OMP_GET_THREAD_NUM()
    !! allocate result array for each thread
    allocate( site_dh(1:fptr% nat), source=0.0 )
    allocate( site_template(1:fptr% nat), source=0 )
    !$OMP DO PRIVATE(nat_loc, typ_loc, coords_loc, idx_loc, tmplt, n, t1, t2, c1, c2, scale, &
    !$OMP            dh_old, dh, t ) SCHEDULE( dynamic, 100)
    do i = 1, fptr% nat
    ! do i = 293855, 293856
       call extract_elements( &
            size(fptr% neighlist, 1), size( fptr% neighlist, 2), fptr% neighlist,&
            size(fptr% veclist, 2), fptr% veclist, fptr% nat, fptr% count_n,&
            i, nat_loc, idx_loc, coords_loc )

       allocate( typ_loc(1:nat_loc) )
       do j = 1, nat_loc
          typ_loc(j) = fptr% typ( idx_loc(j))
       end do

       if( i .eq. 264298 ) then
          write(*,*) nat_loc
          write(*,*) i
          do j = 1, nat_loc
             write(*,*) typ_loc(j), coords_loc(:,j)
          end do
       end if


       ! write(*,*) "on i:", i

       dh_old = 99.9
       t = 0
       do j = 1, fptr% ntemplate
          nullify( tmplt )
          tmplt => fptr% get_template( j )

          if( tmplt% nat .ne. nat_loc ) then
             ! write(*,*) j, tmplt% nat, nat_loc, "skipping"
             cycle
          end if


          n = nat_loc

          allocate( t1, source=tmplt% typ)
          allocate( t2(1:n), source=typ_loc(1:n) )
          if( tmplt% ignore_chem ) then
             !! overwrite atomic types
             t1(:) = 1; t2(:) = 1
          end if

          allocate( c1, source=tmplt% coords )
          allocate( c2, source = coords_loc(:,1:n))
          if( tmplt% rescale ) then
             !! find scale of current conf
             scale = struc_get_scale( n, c2 )
             !! rescale the template
             call struc_rescale( tmplt% nat, c1, scale )
          end if

          dh = struc_get_dh( &
               tmplt% nat, t1, c1, &
               n, t2, c2, &
               kmax_factor, dthr, .true. )

          ! if(me.eq.0) write(*,*) "test", i, j, dh, fptr% ntemplate
          if( i .eq. 264298)write(*,*) "try",i,j,dh

          if( dh .lt. dh_old ) then
             dh_old = dh
             t = j
          end if

          deallocate( t1, t2 )
          deallocate( c1, c2 )

          if( dh_old .le. dthr ) exit
       end do

       if( dh_old < dthr ) then
          !! !$OMP CRITICAL
          fptr% site_dh(i) = dh_old
          fptr% site_template(i) = t
          fptr% site_pg(i) = tmplt% pg
          !! !$OMP END CRITICAL
       end if

       if( i .eq. 264298 ) write(*,*) "HERE",i,dh_old,t

       ! if( me.eq.1)write(*,*) i, dh_old, t
       if( fptr% site_template(i) .eq.0) then
          !$OMP CRITICAL
          block
            type( c_ptr ) :: ccoords, ctyp
            type( c_ptr ) :: cmode
            real( c_double ), pointer :: r2d(:,:)
            integer( c_int ), pointer :: i1d(:)
            integer( c_int ) :: cerr
            allocate( r2d, source=coords_loc)
            ccoords=c_loc(r2d(1,1))
            allocate( i1d, source=typ_loc )
            ctyp = c_loc( i1d(1))

            cmode = f2c_string("nn")
            cerr = citm_add_template( cptr, int(nat_loc, c_int), ctyp, ccoords, &
                 logical(.false., c_bool), cmode, logical(.false., c_bool) )
            deallocate( i1d, r2d )
          end block
          ! write(*,*) "added",fptr% ntemplate, me
          fptr% site_template(i) = fptr% ntemplate
          ! fptr% site_dh = 0.0
          !$OMP END CRITICAL
       end if
       ! if( fptr% site_template(i) .eq. 0 ) write(*,*) i, "still 0"

       ! write(*,*) i, t, dh_old
       ! matched_template = match_site( fptr% template_list, i, dthr, dh )

       deallocate( typ_loc, coords_loc, idx_loc )
    end do
    !$OMP END DO

    deallocate( site_dh, site_template )
    !$OMP END PARALLEL

    ! fptr% site_dh(1)=1.23
    ! fptr% site_template(1)=16

  end function itm_compute



  function itm_compute2( cptr, cthr )result( cerr )bind(C, name="itm_compute2" )
    use m_itm_compute
    implicit none
    type( c_ptr ), value :: cptr
    real( c_double ), value :: cthr
    integer( c_int ) :: cerr

    type( t_itm_ptr ), pointer :: fptr
    real :: dthr


    call c_f_pointer( cptr, fptr )
    dthr = real( cthr )

    !! initial values for result arrays
    fptr% site_dh = 99.9
    fptr% site_template = -1
    fptr% site_pg = -1

    !!
    !! loop over any existing templates
    !!
    call check_pre_templates( fptr, dthr )


    !!
    !! create list of "rough" templates, aka one-shot cshda
    !!
    call create_rough_list( fptr, dthr )


    !!
    !! compare rough templates to create the list of canonical templates
    !!
    write(*,*) "create canon list"
    call create_canon_list( fptr, dthr )
    write(*,*) "canon list done"


    !!
    !! assign canon templates to sites
    !!
    write(*,*) "call assign canon"
    call assign_canon2site( fptr, dthr )



    cerr = 0_c_int
    write(*,*) "finished"
    ! do i = 1, fptr% ntemplate
    !    write(*,*) i, count( fptr% site_template .eq. i )
    ! end do
  end function itm_compute2

  function itm_get_result( cptr, name, cval ) result( cerr ) bind(C, name="itm_get_result" )
    implicit none
    type( c_ptr ), value :: cptr
    type( c_ptr ), value :: name
    type( c_ptr ), intent(out) :: cval
    integer( c_int ) :: cerr

    type( t_itm_ptr ), pointer :: fptr
    character(:), allocatable :: fname
    integer :: dtyp, drank, dsize, dsize2
    integer(c_int), pointer :: i1d(:)
    real( c_double ), pointer :: r1d(:)
    character(:), allocatable :: onestr
    integer :: i

    call c_f_pointer( cptr, fptr )

    ! write(*,*) "entering get_result"
    allocate( fname, source=c2f_string(name) )
    cerr = -1_c_int

    ! write(*,*) "got name",fname
    ! !! get datatype
    ! dtyp = int( itm_get_dtype( name ) )
    ! if( dtyp == ITM_DTYPE_UNKNOWN ) then
    !    write(*,*) "ITM ERROR: unknown data: ",fname
    !    return
    ! end if

    ! !! get datarank
    ! drank = int( itm_get_drank( name ))

    select case( fname )
    case( "site_template" )
       ! write(*,*) "here a"
       allocate( i1d, source=fptr% site_template )
       ! write(*,*) i1d(1)
       cval = c_loc( i1d(1) )
    case( "site_dh" )
       ! write(*,*) "here b"
       allocate( r1d, source=fptr% site_dh )
       ! write(*,*) r1d(1)
       cval = c_loc( r1d(1) )
    case( "site_pg" )
       ! allocate( ipg(1:fptr% nat, source=0))
       ! do i = 1, fptr% nat
       !    ipg(i) = pg_char2int( fptr% site_pg(i) )
       ! end do
       allocate( i1d, source=fptr% site_pg )
       cval = c_loc( i1d(1) )

    case default
       write(*,*) "unknown name:", fname
       return

    end select

    deallocate( fname )
    cerr = 0_c_int
  end function itm_get_result


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

    real :: scale, candidate_match_dh, dh, fthr
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
              kmax_factor, thr, .true. )

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


  subroutine itm_print(cptr, cidx )bind(C,name="itm_print")
    !! print all templates in the list
    implicit none
    type( c_ptr ), value :: cptr
    integer( c_int ), value :: cidx
    type( t_itm_ptr ), pointer :: fptr

    type( t_template ), pointer :: t
    integer :: i, idx

    call c_f_pointer( cptr, fptr )

    idx = int( cidx )

    do i = 1, fptr% ntemplate
       if( i == idx .or. idx == -1 ) then
          t => fptr% get_template( i )
          write(*,*) "got node",i
          call t% print()
       end if
    end do

  end subroutine itm_print


  subroutine itm_print_canon(cptr, cidx )bind(C,name="itm_print_canon")
    !! print all templates in the list
    implicit none
    type( c_ptr ), value :: cptr
    integer( c_int ), value :: cidx
    type( t_itm_ptr ), pointer :: fptr

    type( t_canon ), pointer :: c
    integer :: i, idx

    call c_f_pointer( cptr, fptr )

    idx = int( cidx )

    do i = 1, fptr% ncanon
       if( i == idx .or. idx == -1 ) then
          c => fptr% get_canon( i )
          write(*,*) "got node",i
          call c% print()
       end if
    end do
  end subroutine itm_print_canon


  subroutine itm_check_fast( cptr )bind(C,name="itm_check_fast")
    !! check among the fast templates and match them
    use m_itm_core, only: struc_get_dh
    use omp_lib
    implicit none
    type( c_ptr ), value :: cptr

    type( t_itm_ptr ), pointer :: fptr
    type( t_template ), pointer :: fast1, fast2
    integer :: i, j, idx_i, idx_j
    real :: dh, dthr
    integer, allocatable :: map(:,:)

    call c_f_pointer( cptr, fptr )

    dthr = fptr% dthr

    allocate( map(1:fptr% ntemplate, 1:fptr% ntemplate), source=0)
    do i = 1, fptr% ntemplate
       map(i,i) = 1
    end do

    !$OMP PARALLEL
    !$OMP DO PRIVATE( fast1, fast2, dh, i,j  )
    do i = 1, fptr% ntemplate
       fast1 => fptr% get_template( i )
       do j = i+1, fptr% ntemplate
          fast2 => fptr% get_template( j )

          ! write(*,*) "now comparing:",i,j

          if( fast1% nat .ne. fast2% nat ) then
             ! write(*,*) "different nat",fast1% nat, fast2% nat
             cycle
          end if

          dh = struc_get_dh( fast1% nat, fast1% typ, fast1% coords, &
               fast2% nat, fast2% typ, fast2% coords, kmax_factor, dthr, .false. )

          write(*,*) i,j,dh
          if( dh <= dthr ) then

             ! write(*,*) i, j, dh
!             !$OMP CRITICAL
             map( i, j ) = 1
             map( j, i ) = 1
!             !$OMP END CRITICAL
          end if

       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL


    !$OMP PARALLEL
    !$OMP DO PRIVATE( idx_i, idx_j, fast1 )
    do i = 1, fptr% nat
       !! fast t of this struc
       idx_i = fptr% site_template(i)
       !! minimal index fast of this class
       idx_j = find_first( fptr% ntemplate, map(:,idx_i) )

       if( i .eq. 264298 ) then
          do j = 1, fptr% ntemplate
             if( map(j, idx_i) .eq. 0 ) cycle
             write(*,'(20(i4,x,:))', advance="no") j
          end do
          write(*,*)

          ! write(*,'(*(I1,x,:))') map(:,idx_i)
          write(*,*) i, idx_i, idx_j
       end if

       fptr% site_template(i) = idx_j

       fast1 => fptr% get_template( idx_j )
       fptr% site_pg(i) = fast1% pg

       !! if idx has changed, recompute dh to new template
       if( idx_i /= idx_j ) then
          !! 
       end if

    end do
    !$OMP END DO
    !$OMP END PARALLEL



    deallocate( map)

    ! write(*,*) fptr% nat
    ! write(*,*)
    ! do i = 1, fptr% nat
    !    write(*,*) i, fptr% site_template(i), fptr% site_pg(i)
    ! end do
    ! write(*,*) "exiting check_fast"

  end subroutine itm_check_fast


  subroutine itm_compare_templates(cptr, cidx1, cidx2 )bind(C,name="itm_compare_templates")
    use m_itm_core, only: struc_match_print
    implicit none
    type( c_ptr ), value :: cptr
    integer( c_int ), value :: cidx1, cidx2

    type( t_itm_ptr ), pointer :: fptr
    integer :: idx1, idx2
    type( t_template ), pointer :: t1, t2

    call c_f_pointer( cptr, fptr )
    idx1 = int( cidx1 ); idx2 = int( cidx2 )

    t1 => fptr% get_template( idx1 )
    t2 => fptr% get_template( idx2 )


    call struc_match_print( t1% nat, t1% typ, t1% coords, &
         t2% nat, t2% typ, t2% coords, kmax_factor, fptr% dthr )
  end subroutine itm_compare_templates


  subroutine itm_compare_sites( cptr, csite1, csite2 )bind(C, name="itm_compare_sites")
    implicit none
    type( c_ptr ), value :: cptr
    integer( c_int ), value :: csite1, csite2

    type( t_itm_ptr ), pointer :: fptr
    integer :: isite1, isite2
    integer :: nat1, nat2
    integer, allocatable :: typ1(:), typ2(:)
    real, allocatable :: coords1(:,:), coords2(:,:)
    integer :: i
    real :: rmat(3,3), tr(3), r(3)
    integer, allocatable :: perm(:)
    real :: dh
    integer, allocatable :: c1(:), c2(:)
    integer :: ierr


    call c_f_pointer( cptr, fptr )

    isite1 = int( csite1 )
    isite2 = int( csite2 )

    call get_local_conf(fptr, isite1, nat1, typ1, coords1 )
    call get_local_conf(fptr, isite2, nat2, typ2, coords2 )

    allocate( perm(1:nat2))



    if( nat1 /= nat2 ) then
       write(*,*) "number of atoms not equal"
       return
    end if

    write(*,*) nat1
    write(*,*) "original isite:",isite1
    do i = 1, nat1
       write(*,*) typ1(i), coords1(:,i)
    end do


    allocate( c1(1:nat1))
    allocate( c2(1:nat2))
    c1 = 0; c1(1) = 1
    c2 = 0; c2(1) = 1

    call ira_unify( nat1, typ1, coords1, c2, &
         nat2, typ2, coords2, c2, 1.8, rmat, tr, perm, dh, ierr )

    deallocate( c1, c2)

    !! permute
    typ2 = typ2(perm)
    coords2(:,:) = coords2(:,perm)

    call svdrot_m( nat1, typ1, coords1, nat2, typ2, coords2, rmat, tr, ierr )

    !! apply
    do i = 1, nat2
       coords2(:,i) = matmul(rmat, coords2(:,i)) + tr
    end do

    write(*,*) nat2
    write(*,*) "matched isite:",isite2
    do i = 1, nat2
       write(*,*) typ2(i), coords2(:,i)
    end do


    write(*,*) "initial ira:",dh
    dh = 0.0
    do i = 1, nat1
       r = coords2(:,i) - coords1(:,i)
       dh = max( dh, norm2(r))
    end do
    write(*,*) "final dh:",dh

    deallocate( typ1, coords1 )
    deallocate( typ2, coords2 )
  end subroutine itm_compare_sites
  subroutine itm_compare_site_template( cptr, csite, ctmplt )bind(C,name="itm_compare_site_template")
    implicit none
    type( c_ptr ), value :: cptr
    integer( c_int ), value :: csite, ctmplt

    type( t_itm_ptr ), pointer :: fptr
    integer :: isite, itmplt
    type( t_template ), pointer :: t
    integer :: nat_loc
    integer, allocatable :: typ_loc(:), typ2(:)
    real, allocatable :: coords_loc(:,:), coords2(:,:)
    real :: rmat(3,3), tr(3), r(3)
    integer, allocatable :: perm(:), c1(:), c2(:)
    real :: dh
    integer :: ierr, i


    isite = int( csite )
    itmplt = int( ctmplt )

    call c_f_pointer( cptr, fptr )

    !! get the template
    t => fptr% get_template( itmplt )

    call get_local_conf( fptr, isite, nat_loc, typ_loc, coords_loc )


    write(*,*) t% nat
    write(*,*) "original template index:", itmplt
    do i = 1, t% nat
       write(*,*) t% typ(i), t% coords(:,i)
    end do



    allocate( perm(1:t% nat))
    dh = dh_cshda( nat_loc, typ_loc, coords_loc, t% nat, t% typ, t% coords, 999.9, perm )
    write(*,*) "initial cshda:", dh

    write(*,*) nat_loc
    write(*,*) "original site index:",isite
    do i = 1, nat_loc
       write(*,*) typ_loc(i), coords_loc(:,i)
    end do

    if( nat_loc /= t% nat) then
       write(*,*) "number of atoms not equal"
       return
    end if


    allocate( c1(1:nat_loc))
    allocate( c2(1:t% nat))
    c1 = 0; c1(1) = 1
    c2 = 0; c2(1) = 1

    call ira_unify( nat_loc, typ_loc, coords_loc, c2, &
         t% nat, t% typ, t% coords, c2, 1.8, rmat, tr, perm, dh, ierr )
    deallocate( c1, c2)

    !! permute
    allocate( typ2, source=t% typ(perm))
    allocate( coords2, source=t% coords(:,perm))

    call svdrot_m( nat_loc, typ_loc, coords_loc, t% nat, typ2, coords2, rmat, tr, ierr )

    !! apply
    do i = 1, t% nat
       coords2(:,i) = matmul(rmat, coords2(:,i)) + tr
    end do

    write(*,*) t% nat
    write(*,*) "matched template"
    do i = 1, t% nat
       write(*,*) typ2(i), coords2(:,i)
    end do


    write(*,*) "initial ira:",dh
    dh = 0.0
    do i = 1, nat_loc
       r = coords2(:,i) - coords_loc(:,i)
       dh = max( dh, norm2(r))
    end do
    write(*,*) "final dh:",dh

    deallocate( typ2, coords2 )
    nullify(t)
  end subroutine itm_compare_site_template


  subroutine itm_compare_site_canon( cptr, csite, ccanon )bind(C,name="itm_compare_site_canon")
    implicit none
    type( c_ptr ), value :: cptr
    integer( c_int ), value :: csite, ccanon

    type( t_itm_ptr ), pointer :: fptr
    integer :: isite, icanon
    type( t_canon ), pointer :: c
    integer :: nat_loc
    integer, allocatable :: typ_loc(:), typ2(:)
    real, allocatable :: coords_loc(:,:), coords2(:,:)
    real :: rmat(3,3), tr(3), r(3)
    integer, allocatable :: perm(:), c1(:), c2(:)
    real :: dh
    integer :: ierr, i


    isite = int( csite )
    icanon = int( ccanon )

    call c_f_pointer( cptr, fptr )

    !! get the canonical template
    c => fptr% get_canon( icanon )

    call get_local_conf( fptr, isite, nat_loc, typ_loc, coords_loc )


    write(*,*) c% nat
    write(*,*) "original canon index:", icanon
    do i = 1, c% nat
       write(*,*) c% typ(i), c% coords(:,i)
    end do



    allocate( perm(1:c% nat))
    dh = dh_cshda( nat_loc, typ_loc, coords_loc, c% nat, c% typ, c% coords, 999.9, perm )
    write(*,*) "initial cshda:", dh

    write(*,*) nat_loc
    write(*,*) "original site index:",isite
    do i = 1, nat_loc
       write(*,*) typ_loc(i), coords_loc(:,i)
    end do

    if( nat_loc /= c% nat) then
       write(*,*) "number of atoms not equal", nat_loc, c% nat
       return
    end if


    allocate( c1(1:nat_loc))
    allocate( c2(1:c% nat))
    c1 = 0; c1(1) = 1
    c2 = 0; c2(1) = 1

    call ira_unify( nat_loc, typ_loc, coords_loc, c2, &
         c% nat, c% typ, c% coords, c2, 1.8, rmat, tr, perm, dh, ierr )
    deallocate( c1, c2)

    !! permute
    allocate( typ2, source=c% typ(perm))
    allocate( coords2, source=c% coords(:,perm))

    call svdrot_m( nat_loc, typ_loc, coords_loc, c% nat, typ2, coords2, rmat, tr, ierr )

    !! apply
    do i = 1, c% nat
       coords2(:,i) = matmul(rmat, coords2(:,i)) + tr
    end do

    write(*,*) c% nat
    write(*,*) "matched canon"
    do i = 1, c% nat
       write(*,*) typ2(i), coords2(:,i)
    end do


    write(*,*) "initial ira:",dh
    dh = 0.0
    do i = 1, nat_loc
       r = coords2(:,i) - coords_loc(:,i)
       dh = max( dh, norm2(r))
    end do
    write(*,*) "final dh:",dh

    deallocate( typ2, coords2 )
    nullify(c)
  end subroutine itm_compare_site_canon


  !!-----------------------------------------------------------
  !! local functions

  function find_first( n, arr )result(idx)
    !! find first nonzero entry in array
    integer, intent(in) :: n
    integer, dimension(n), intent(in) :: arr
    integer :: idx

    integer :: i

    do i = 1, n
       if( arr(i) .eq. 0 ) cycle
       exit
    end do
    idx = i
  end function find_first



  ! function itm_get_dtype( cptr, name )result(ctype)bind(C,name="itm_get_dtype")
  function itm_get_dtype( name )result(ctype)bind(C,name="itm_get_dtype")
    implicit none
    ! type( c_ptr ), value :: cptr
    type( c_ptr ), value :: name
    integer( c_int ) :: ctype

    character(:), allocatable :: fname
    integer :: dtype

    allocate( fname, source=c2f_string(name) )

    dtype = ITM_DTYPE_UNKNOWN
    select case( fname )
    case( "nat", "typ", "neighlist", "count_n", "site_template", "ntemplate", "site_pg" )
       dtype = ITM_DTYPE_INT

    case( "veclist", "site_dh" )
       dtype = ITM_DTYPE_REAL

    case default
       dtype = -1
    end select

    deallocate( fname )
    ctype = int( dtype, c_int )
  end function itm_get_dtype

  ! function itm_get_drank( cptr, name )result(crank)bind(C,name="itm_get_drank")
  function itm_get_drank( name )result(crank)bind(C,name="itm_get_drank")
    implicit none
    ! type( c_ptr ), value :: cptr
    type( c_ptr ), value :: name
    integer( c_int ) :: crank

    character(:), allocatable :: fname
    integer :: drank

    allocate( fname, source=c2f_string(name) )

    select case( fname )
    case( "nat", "ntemplate" )
       drank = 0

    case( "site_template", "site_dh", "typ", "count_n", "site_pg" )
       drank = 1

    case( "neighlist", "veclist" )
       drank = 2

    case default
       drank = -1
    end select

    deallocate( fname )
    crank = int( drank, c_int )
  end function itm_get_drank

  function itm_get_dsize( cptr, name )result(csize)bind(C,name="itm_get_dsize")
    !! useful when extracting data
    implicit none
    type( c_ptr ), value :: cptr
    type( c_ptr ), value :: name
    type( c_ptr ) :: csize

    type( t_itm_ptr ), pointer :: fptr
    integer( c_int ), pointer :: psize(:)
    integer, allocatable :: fsize(:)

    character(:), allocatable :: fname

    call c_f_pointer( cptr, fptr )

    allocate( fname, source=c2f_string(name) )

    select case( fname )
    case( "site_template", "site_dh", "site_pg" )
       allocate( fsize(1), source = fptr% nat )
       allocate( psize, source=fsize )
       csize = c_loc( psize(1) )
       deallocate( fsize )
    end select

    deallocate( fname )
  end function itm_get_dsize

  function itm_pg_int2char( pgint )result(pgchar)bind(C,name="itm_pg_int2char")
    use m_itm_core, only: pg_int2char
    implicit none
    integer( c_int ), value :: pgint
    type( c_ptr ) :: pgchar

    character(len=5) :: pg

    pg = pg_int2char( pgint )
    pgchar = f2c_string(pg)

  end function itm_pg_int2char


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
  FUNCTION f2c_string(f_string) RESULT(ptr)
    use iso_c_binding
    interface
       FUNCTION c_malloc(size) BIND(C, name='malloc')
         IMPORT :: c_ptr, c_size_t
         IMPLICIT NONE
         INTEGER(c_size_t), VALUE :: size
         TYPE(c_ptr) :: c_malloc
       END FUNCTION c_malloc
    end interface

    CHARACTER(LEN=*), INTENT(IN)           :: f_string
    CHARACTER(LEN=1, KIND=c_char), POINTER :: c_string(:)
    TYPE(c_ptr) :: ptr
    INTEGER(c_size_t) :: i, n

    n = LEN_TRIM(f_string)
    ptr = c_malloc(n+1)
    CALL C_F_POINTER(ptr, c_string, [n+1])
    DO i=1, n
       c_string(i) = f_string(i:i)
    END DO
    c_string(n+1) = c_null_char
  END FUNCTION f2c_string

end module itm
