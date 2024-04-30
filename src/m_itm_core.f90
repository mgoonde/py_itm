module m_itm_core

  !! functions to match, resize, etc

  implicit none


contains

  function struc_get_scale( nat, coords_in )result(scale)
    !! return the scale, as distance from center
    implicit none
    integer, intent(in) :: nat
    real, dimension(3,nat), intent(in) :: coords_in
    real :: scale

    integer :: i, imin
    real :: gc(3), rmin, r, rdum(3)
    real, dimension(3,nat) :: coords

    !! working copy
    coords = coords_in

    !! recenter to geometrical center
    gc = sum( coords(:,:), 2)/nat
    do i = 1, nat
       coords(:,i) = coords(:,i) - gc(:)
    end do

    rmin = 999.9
    !! ensure first atom is nearest to gc
    !! find atom closest to gc, put that atom at first place
    do i = 1, nat
       r = norm2( coords(:,i) )
       if( r .lt. rmin ) then
          rmin = r
          imin = i
       end if
    end do
    !! tmp copy, then switch
    rdum = coords(:,imin)
    coords(:,imin) = coords(:,1)
    coords(:,1) = rdum


    !! get scale
    scale = 0.0
    do i = 1, nat
       scale = scale + norm2( coords(:,i) )
       ! write(*,*) i, coords(:,i), norm2(coords(:,i))
    end do
    scale = scale / nat

  end function struc_get_scale

  subroutine struc_rescale( nat, coords, scale )
    !! rescale coords to scale
    implicit none
    integer, intent(in) :: nat
    real, dimension(3,nat), intent(inout) :: coords
    real, intent(in) :: scale

    integer :: i

    do i = 1, nat
       coords(:,i) = coords(:,i)*scale
    end do

  end subroutine struc_rescale


  function struc_get_dh( nat1, typ1, coords1, nat2, typ2, coords2, kmax_factor, dthr, fast )result(dh)
    implicit none
    integer, intent(in) :: nat1
    integer, dimension(nat1), intent(in) :: typ1
    real, dimension(3,nat1), intent(in) :: coords1
    integer, intent(in) :: nat2
    integer, dimension(nat2), intent(in) :: typ2
    real, dimension(3,nat2), intent(in) :: coords2
    real, intent(in) :: kmax_factor
    real, intent(in) :: dthr
    logical, intent(in) :: fast
    real :: dh, rmsd

    integer, dimension(nat1) :: c1
    integer, dimension(nat2) :: c2
    real, dimension(3) :: tr, str, rdum
    real, dimension(3,3) :: rmat, srmat
    integer, dimension(nat1) :: p
    integer :: i, ierr, serr
    character(512) :: msg

    integer, dimension(nat2) :: found
    real, dimension(nat2) :: dists
    real, dimension(3,nat1) :: flipx
    integer, dimension(nat1) :: typ1_w
    real, dimension(3,nat1) :: coords1_w
    integer, dimension(nat2) :: typ2_w
    real, dimension(3,nat2) :: coords2_w


#ifdef DEBUG
    write(*,*) "struc_get_dh:: now matching:"
    write(*,*) nat1
    write(*,*) '1'
    do i = 1, nat1
       write(*,*) typ1(i), coords1(:,i)
    end do
    write(*,*) nat2
    write(*,*) '2'
    do i = 1, nat2
       write(*,*) typ2(i), coords2(:,i)
    end do
#endif



    !! work arrays
    typ1_w = typ1; coords1_w = coords1
    typ2_w = typ2; coords2_w = coords2


    dh = 999.9
    ! ierr = 0

    !! check on number fo atoms, this should never happen
    ! if( nat1 /= nat2 ) then
    !    write(*,*) "in struc_get_dh, nat not equal!",nat1,nat2
    !    return
    ! end if

    !! set candidates to first index
    c1(:) = 0; c1(1) = 1
    c2(:) = 0; c2(1) = 1


    !! first try on cshda
    call cshda( nat1, typ1_w, coords1_w, nat2, typ2_w, coords2_w, 1.1*dthr, found, dists )
    dh = maxval( dists )
    if( dh .lt. dthr ) return


    !! fast mode
    if( fast ) return


    !! flip the x, try again (man, this is arbitrary...!)
    do i = 1, nat1
       coords1_w(1,i) = -coords1_w(1,i)
    end do
    call cshda( nat1, typ1_w, coords1_w, nat2, typ2_w, coords2_w, 1.1*dthr, found, dists )
    dh = maxval( dists )
    if( dh .lt. dthr ) return
    !! reset
    coords1_w = coords1

    !! get first svd

    !! recompute cshda


    !! call ira to obtain dh
#ifdef DEBUG
    write(*,*) "full ira"
#endif DEBUG

    call ira_unify( nat1, typ1_w, coords1_w, c1, &
                    nat2, typ2_w, coords2_w, c2, &
                    kmax_factor, rmat, tr, p, dh, ierr )
    if( ierr /= 0 ) then
       write(*,*) "err"
       call ira_get_errmsg( ierr, msg )
       write(*,*) trim(msg)
    end if

#ifdef DEBUG
    write(*,*) dh, rmsd
    write(*,*) "rmat from ira_unify"
    write(*,'(3f8.4)') rmat(1,:)
    write(*,'(3f8.4)') rmat(2,:)
    write(*,'(3f8.4)') rmat(3,:)
    write(*,*) "tr from ira_unify"
    write(*,'(3f8.4)') tr
    write(*,*)
#endif



    typ2_w = typ2_w(p)
    coords2_w = coords2_w(:,p)
    do i = 1, nat2
       coords2_w(:,i) = matmul( rmat, coords2_w(:,i) ) + tr
    end do
    !! call svd
    call svdrot_m( nat1, typ1_w, coords1_w, &
         nat2, typ2_w(1:nat1), coords2_w(:,1:nat1), srmat, str, ierr )
    if( ierr /= 0 ) then
       write(*,*) "err svd"
       return
    end if

#ifdef DEBUG
    write(*,*) "rmat from svd"
    write(*,'(3f8.4)') srmat(1,:)
    write(*,'(3f8.4)') srmat(2,:)
    write(*,'(3f8.4)') srmat(3,:)
    write(*,*) "tr from svd"
    write(*,'(3f8.4)') str
    write(*,*)
#endif


    rmat = matmul( srmat, rmat )
    tr = matmul(srmat, tr) + str
    ! reset
    typ2_w = typ2(p)
    coords2_w = coords2(:,p)
    ! apply svd
    do i = 1, nat2
       coords2_w(:,i) = matmul( rmat, coords2_w(:,i) ) + tr
    end do

    ! get dh
    dh = 0.0
    do i = 1, nat1
       rdum = coords2_w(:,i) - coords1(:,i)
       dh = max( dh, sqrt(dot_product(rdum,rdum)))
       ! dh = dh + dot_product( rdum, rdum )
    end do
    ! dh = sqrt(dh)/nat1


#ifdef DEBUG
    write(*,*) nat1+nat2
    write(*,*) 'matched'
    do i = 1, nat1
       write(*,*) 1, coords1(:,i)
       write(*,*) 2, coords2_w(:,i)
    end do
#endif




  end function struc_get_dh


  subroutine struc_match_print( nat1, typ1, coords1, nat2, typ2, coords2, kmax_factor, dthr )
    implicit none
    integer, intent(in) :: nat1
    integer, dimension(nat1), intent(in) :: typ1
    real, dimension(3,nat1), intent(in) :: coords1
    integer, intent(in) :: nat2
    integer, dimension(nat2), intent(in) :: typ2
    real, dimension(3,nat2), intent(in) :: coords2
    real, intent(in) :: kmax_factor
    real, intent(in) :: dthr
    real :: dh, rmsd

    integer, dimension(nat1) :: c1
    integer, dimension(nat2) :: c2
    real, dimension(3) :: tr, str, rdum
    real, dimension(3,3) :: rmat, srmat
    integer, dimension(nat1) :: p
    integer :: i, ierr, serr
    character(512) :: msg

    integer, dimension(nat2) :: found
    real, dimension(nat2) :: dists
    real, dimension(3,nat1) :: flipx
    integer, dimension(nat1) :: typ1_w
    real, dimension(3,nat1) :: coords1_w
    integer, dimension(nat2) :: typ2_w
    real, dimension(3,nat2) :: coords2_w


    !! work arrays
    typ1_w = typ1; coords1_w = coords1
    typ2_w = typ2; coords2_w = coords2


    dh = 999.9

    !! set candidates to first index
    c1(:) = 0; c1(1) = 1
    c2(:) = 0; c2(1) = 1

    !! error for cshda
    if( nat1 /= nat2 ) return

    call cshda( nat1, typ1_w, coords1_w, nat2, typ2_w, coords2_w, 999.9, found, dists )
    write(*,*) "initial csdha dh:",maxval(dists)

    call ira_unify( nat1, typ1_w, coords1_w, c1, &
                    nat2, typ2_w, coords2_w, c2, &
                    kmax_factor, rmat, tr, p, dh, ierr )
    if( ierr /= 0 ) then
       write(*,*) "err"
       call ira_get_errmsg( ierr, msg )
       write(*,*) trim(msg)
    end if

    typ2_w = typ2_w(p)
    coords2_w = coords2_w(:,p)
    do i = 1, nat2
       coords2_w(:,i) = matmul( rmat, coords2_w(:,i) ) + tr
    end do
    !! call svd
    call svdrot_m( nat1, typ1_w, coords1_w, &
         nat2, typ2_w(1:nat1), coords2_w(:,1:nat1), srmat, str, ierr )
    if( ierr /= 0 ) then
       write(*,*) "err svd"
       return
    end if

#ifdef DEBUG
    write(*,*) "rmat from svd"
    write(*,'(3f8.4)') srmat(1,:)
    write(*,'(3f8.4)') srmat(2,:)
    write(*,'(3f8.4)') srmat(3,:)
    write(*,*) "tr from svd"
    write(*,'(3f8.4)') str
    write(*,*)
#endif


    rmat = matmul( srmat, rmat )
    tr = matmul(srmat, tr) + str
    ! reset
    typ2_w = typ2(p)
    coords2_w = coords2(:,p)
    ! apply svd
    do i = 1, nat2
       coords2_w(:,i) = matmul( rmat, coords2_w(:,i) ) + tr
    end do

    ! get dh
    dh = 0.0
    do i = 1, nat1
       rdum = coords2_w(:,i) - coords1(:,i)
       dh = max( dh, sqrt(dot_product(rdum,rdum)))
       ! dh = dh + dot_product( rdum, rdum )
    end do
    ! dh = sqrt(dh)/nat1


    write(*,*) nat1
    write(*,*) 't1, matched dh=',dh
    do i = 1, nat1
       write(*,*) typ1(i), coords1(:,i)
    end do
    write(*,*) nat2
    write(*,*) "t2"
    do i = 1, nat2
       write(*,*) typ2_w(i), coords2_w(:,i)
    end do
  end subroutine struc_match_print



  subroutine order_atoms( nat, typ, coords )
    !! order atoms by distance from origin
    implicit none
    integer, intent(in) :: nat
    integer, dimension(nat), intent(inout) :: typ
    real, dimension(3,nat), intent(inout) :: coords

    real, dimension(2,nat) :: d_o  !! distance, order
    integer :: i

    do i = 1, nat
       d_o(1,i) = norm2( coords(:,i) )
       d_o(2,i) = real(i)
    end do

    !! sort d_o by dim 1:
    !! routine in IRA ibrary
    call sort( nat, 2, d_o, 1 )

    !! permute coords by d_o dim 2
    coords(:,:) = coords(:, nint(d_o(2,:)) )
    typ(:) = typ(nint(d_o(2,:)) )

  end subroutine order_atoms


  function find_dmax( nat, coords )result( dmax )
    !! find the largest distance from the first atom
    implicit none
    integer, intent(in) :: nat
    real, dimension(3,nat), intent(in) :: coords
    real :: dmax

    integer :: i
    real, dimension(3) :: rdum

    dmax = 0.0
    do i = 1, nat
       rdum = coords(:,1) - coords(:,i)
       dmax = max( dmax, sqrt(dot_product(rdum, rdum)) )
    end do
    return
  end function find_dmax

  function natoms_dist_less_than( nat, coords, rcut ) result( n )
    !! return number of atoms whose distance from origin is less than rcut.
    !! The coords have to be ordered by distance fromthe origin on input!
    implicit none
    integer, intent(in) :: nat
    real, dimension(3,nat), intent(in) :: coords
    real, intent(in) :: rcut
    integer :: n

    integer :: i

    !! error rcut value
    n = -1
    if( rcut  .lt. 1e-8) return

    n = 0
    do i = 1, nat
       if( norm2(coords(:,i)) .gt. rcut) exit
       n = n + 1
    end do
  end function natoms_dist_less_than


  function count_unique( n, a_in )result( n_unique )
    !! count unique values (larger than -99) in integer 1d array
    implicit none
    integer, intent(in) :: n
    integer, dimension(n), intent(in) :: a_in
    integer :: n_unique

    integer, dimension(n) :: a
    integer :: i, j, k

    n_unique = 0

    !! working copy
    a = a_in

    !! first value in a
    k = a(1)
    do i = 1, n
       if( k .eq. -99 ) exit
       !! add count of unique
       n_unique = n_unique + 1

       !! put all values where a==k to -99
       do j = 1, n
          if( a(j) == k) a(j) = -99
       end do

       !! set next k
       k = maxval( a, 1)
    end do


  end function count_unique


  subroutine count_nn( dim1, dim2, neighlist, nat, count_n )
    !! get the number of neighbours for each atom,
    !! the neighlist has entries like: neighlist(j) = ( 3, 462 )
    !! where the first index is central atom, and second is the neighbor index.
    !! This routine first counts how many neighbors each atom has, then it
    !! makes a cumulative sum, which makes it faster to parse the neigharray later on.
    !!
    !! Using cumulative sum, the value count_n( target_idx-1 ) is the starting index in the
    !! neighlist array, and the value count_n( target_idx ) is the final index for atom index target_idx
    implicit none
    integer, intent(in) :: dim1
    integer, intent(in) :: dim2
    integer, dimension(dim1, dim2), intent(in) :: neighlist
    integer, intent(in) :: nat
    integer, dimension(nat), intent(out) :: count_n

    integer :: me, count_me, i, j

    !! neighlist is sorted on the first index, count how many entires
    !! have the same first index (this is number of neighbors for atom idx1)
    me = 1
    count_me = 0
    do i = 1, dim2
       !! first index is me
       j = neighlist(1,i)
       !! first index changed, write data and reset counter
       if( j .ne. me ) then
          count_n(me) = count_me
          me = j
          count_me = 0
       end if
       !! else increase counter
       count_me = count_me + 1
    end do
    !! the last entry
    count_n(me) = count_me

    !! now make this a cumulative sum
    do i = 2, nat
       count_n(i) = count_n(i-1) + count_n(i)
    end do

  end subroutine count_nn


  subroutine extract_elements( dim1, dim2, neighlist, vldim, veclist, nat, count_n, tgt_idx, nloc, idx, vec )
    !! extract neighbours and vectors for atom index tgt_idx
    implicit none
    integer, intent(in) :: dim1, dim2
    integer, dimension(dim1, dim2), intent(in) :: neighlist
    integer, intent(in) :: vldim
    real, dimension(3,vldim), intent(in) :: veclist
    integer, intent(in) :: nat
    integer, dimension(nat), intent(in) :: count_n
    integer, intent(in) :: tgt_idx
    integer, intent(out) :: nloc
    integer, allocatable, intent(out) :: idx(:)
    real, allocatable, intent(out) :: vec(:,:)

    integer :: si, ei, i, j

    !! end index
    ei = count_n( tgt_idx )

    !! start index
    if( tgt_idx == 1 ) then
       si = 1
       !! number of atoms in local conf: number of neighbours plus one for central
       nloc = ei + 1
    else
       si = count_n( tgt_idx - 1 ) + 1
       nloc = ei - si + 2
    end if

    ! write(*,*) "starting:",si
    ! write(*,*) "ending:",ei
    ! write(*,*) "nloc",nloc

    ! write(*,*) "got"
    ! do i = si, ei
    !    write(*,*) i
    ! end do
    ! write(*,*) '='

    allocate( idx(1:nloc) )
    allocate( vec(1:3,1:nloc))

    !! first vector is zero
    vec(:,1) = (/0.0, 0.0, 0.0/)
    j = 1
    ! write(*,*) size(vec(:,2:),2)
    ! write(*,*) size( veclist(:,si:ei),2)
    vec(:,2:) = veclist(:,si:ei)
    ! do i = si, ei
    !    ! write(*,*) i, veclist(:,i)
    !    j = j + 1
    !    vec(:,j) = veclist(:, i )
    ! end do

    idx(1) = tgt_idx
    idx(2:) = neighlist(2, si:ei)
  end subroutine extract_elements


  function pg_char2int( pgchar )result(pgint)
    !! convert pg string into integer encoder
    implicit none
    character(*), intent(in) :: pgchar
    integer :: pgint

    integer :: i
    integer, parameter :: npg=204
    character(len=5), dimension(npg) :: pgdata = [&
         '   C1', '   Cs', '   C2', '  C2h', '  C2v', '   C3', '  C3h', '  C3v', '   C4', '  C4h', '  C4v', &
         '   C5', '  C5h', '  C5v', '   C6', '  C6h', '  C6v', '   C7', '  C7h', '  C7v', '   C8', '  C8h', &
         '  C8v', '   C9', '  C9h', '  C9v', '  C10', ' C10h', ' C10v', '   D2', '  D2h', '  D2d', '   D3', &
         '  D3h', '  D3d', '   D4', '  D4h', '  D4d', '   D5', '  D5h', '  D5d', '   D6', '  D6h', '  D6d', &
         '   D7', '  D7h', '  D7d', '   D8', '  D8h', '  D8d', '   D9', '  D9h', '  D9d', '  D10', ' D10h', &
         ' D10d', '    T', '   Th', '   Td', '    O', '   Oh', '    I', '   Ih', '   Ci', '   S4', '   S6', &
         '   S8', '  S10', &
         '  C1+', '  Cs+', '  C2+', ' C2h+', ' C2v+', '  C3+', ' C3h+', ' C3v+', '  C4+', ' C4h+', ' C4v+', &
         '  C5+', ' C5h+', ' C5v+', '  C6+', ' C6h+', ' C6v+', '  C7+', ' C7h+', ' C7v+', '  C8+', ' C8h+', &
         ' C8v+', '  C9+', ' C9h+', ' C9v+', ' C10+', 'C10h+', 'C10v+', '  D2+', ' D2h+', ' D2d+', '  D3+', &
         ' D3h+', ' D3d+', '  D4+', ' D4h+', ' D4d+', '  D5+', ' D5h+', ' D5d+', '  D6+', ' D6h+', ' D6d+', &
         '  D7+', ' D7h+', ' D7d+', '  D8+', ' D8h+', ' D8d+', '  D9+', ' D9h+', ' D9d+', ' D10+', 'D10h+', &
         'D10d+', '   T+', '  Th+', '  Td+', '   O+', '  Oh+', '   I+', '  Ih+', '  Ci+', '  S4+', '  S6+', &
         '  S8+', ' S10+', &
         '  C1-', '  Cs-', '  C2-', ' C2h-', ' C2v-', '  C3-', ' C3h-', ' C3v-', '  C4-', ' C4h-', ' C4v-', &
         '  C5-', ' C5h-', ' C5v-', '  C6-', ' C6h-', ' C6v-', '  C7-', ' C7h-', ' C7v-', '  C8-', ' C8h-', &
         ' C8v-', '  C9-', ' C9h-', ' C9v-', ' C10-', 'C10h-', 'C10v-', '  D2-', ' D2h-', ' D2d-', '  D3-', &
         ' D3h-', ' D3d-', '  D4-', ' D4h-', ' D4d-', '  D5-', ' D5h-', ' D5d-', '  D6-', ' D6h-', ' D6d-', &
         '  D7-', ' D7h-', ' D7d-', '  D8-', ' D8h-', ' D8d-', '  D9-', ' D9h-', ' D9d-', ' D10-', 'D10h-', &
         'D10d-', '   T-', '  Th-', '  Td-', '   O-', '  Oh-', '   I-', '  Ih-', '  Ci-', '  S4-', '  S6-', &
         '  S8-', ' S10-' ]

    do i = 1, npg
       if( trim(adjustl(pgchar)) == trim(adjustl(pgdata(i))) ) then
          pgint = i
          exit
       end if
    end do
  end function pg_char2int

  function pg_int2char( pgint )result( pgchar )
    !! convert pg integer encoder into string
    implicit none
    integer, intent(in) :: pgint
    character(len=5) :: pgchar
    integer, parameter :: npg=204
    character(len=5), dimension(npg) :: pgdata = [&
         '   C1', '   Cs', '   C2', '  C2h', '  C2v', '   C3', '  C3h', '  C3v', '   C4', '  C4h', '  C4v', &
         '   C5', '  C5h', '  C5v', '   C6', '  C6h', '  C6v', '   C7', '  C7h', '  C7v', '   C8', '  C8h', &
         '  C8v', '   C9', '  C9h', '  C9v', '  C10', ' C10h', ' C10v', '   D2', '  D2h', '  D2d', '   D3', &
         '  D3h', '  D3d', '   D4', '  D4h', '  D4d', '   D5', '  D5h', '  D5d', '   D6', '  D6h', '  D6d', &
         '   D7', '  D7h', '  D7d', '   D8', '  D8h', '  D8d', '   D9', '  D9h', '  D9d', '  D10', ' D10h', &
         ' D10d', '    T', '   Th', '   Td', '    O', '   Oh', '    I', '   Ih', '   Ci', '   S4', '   S6', &
         '   S8', '  S10', &
         '  C1+', '  Cs+', '  C2+', ' C2h+', ' C2v+', '  C3+', ' C3h+', ' C3v+', '  C4+', ' C4h+', ' C4v+', &
         '  C5+', ' C5h+', ' C5v+', '  C6+', ' C6h+', ' C6v+', '  C7+', ' C7h+', ' C7v+', '  C8+', ' C8h+', &
         ' C8v+', '  C9+', ' C9h+', ' C9v+', ' C10+', 'C10h+', 'C10v+', '  D2+', ' D2h+', ' D2d+', '  D3+', &
         ' D3h+', ' D3d+', '  D4+', ' D4h+', ' D4d+', '  D5+', ' D5h+', ' D5d+', '  D6+', ' D6h+', ' D6d+', &
         '  D7+', ' D7h+', ' D7d+', '  D8+', ' D8h+', ' D8d+', '  D9+', ' D9h+', ' D9d+', ' D10+', 'D10h+', &
         'D10d+', '   T+', '  Th+', '  Td+', '   O+', '  Oh+', '   I+', '  Ih+', '  Ci+', '  S4+', '  S6+', &
         '  S8+', ' S10+', &
         '  C1-', '  Cs-', '  C2-', ' C2h-', ' C2v-', '  C3-', ' C3h-', ' C3v-', '  C4-', ' C4h-', ' C4v-', &
         '  C5-', ' C5h-', ' C5v-', '  C6-', ' C6h-', ' C6v-', '  C7-', ' C7h-', ' C7v-', '  C8-', ' C8h-', &
         ' C8v-', '  C9-', ' C9h-', ' C9v-', ' C10-', 'C10h-', 'C10v-', '  D2-', ' D2h-', ' D2d-', '  D3-', &
         ' D3h-', ' D3d-', '  D4-', ' D4h-', ' D4d-', '  D5-', ' D5h-', ' D5d-', '  D6-', ' D6h-', ' D6d-', &
         '  D7-', ' D7h-', ' D7d-', '  D8-', ' D8h-', ' D8d-', '  D9-', ' D9h-', ' D9d-', ' D10-', 'D10h-', &
         'D10d-', '   T-', '  Th-', '  Td-', '   O-', '  Oh-', '   I-', '  Ih-', '  Ci-', '  S4-', '  S6-', &
         '  S8-', ' S10-' ]

    pgchar = "none"
    if( pgint > npg .or. pgint < 1 ) return

    pgchar = trim(pgdata(pgint))
  end function pg_int2char


  !> @details
  !! obtain directly dh of two structures from cshda. If not below dthr, dh can be incorrect (huge value)
  function dh_cshda( nat1, typ1, coords1, nat2, typ2, coords2, dthr, perm ) result(dh)
    integer, intent(in) :: nat1
    integer, intent(in) :: typ1( nat1 )
    real,    intent(in) :: coords1( 3, nat1 )
    integer, intent(in) :: nat2
    integer, intent(in) :: typ2( nat2 )
    real,    intent(in) :: coords2( 3, nat2 )
    real,    intent(in) :: dthr
    integer, intent(out) :: perm(nat2)
    real :: dh

    real, allocatable :: dists(:)

    dh = 999.9
    if( nat1 /= nat2 ) return

    allocate( dists(1:nat2))

    call cshda( nat1, typ1, coords1, nat2, typ2, coords2, dthr, perm, dists )

    dh = maxval( dists(1:nat1), 1)
    deallocate( dists )
  end function dh_cshda


  function dh_full_ira( nat1, typ1, coords1, nat2, typ2_in, coords2_in, kmax, perm ) result(dh)
    implicit none
    integer, intent(in) :: nat1
    integer, intent(in) :: typ1(nat1)
    real,    intent(in) :: coords1(3,nat1)
    integer, intent(in) :: nat2
    integer, intent(in) :: typ2_in(nat2)
    real,    intent(in) :: coords2_in(3,nat2)
    real,    intent(in) :: kmax
    integer, intent(out) :: perm(nat2)
    real :: dh

    !! local arrays
    integer :: typ2(nat2)
    real :: coords2(3,nat2)
    integer :: c1(nat1), c2(nat2)
    real :: rmat(3,3), tr(3)
    real :: dx, dy, dz
    integer :: ierr
    character(512) :: msg
    integer :: i

    ! write(*,*) "enter full ira"
    dh = 999.9
    if( nat1 /= nat2 ) return

    !! copy input struc2
    typ2 = typ2_in
    coords2 = coords2_in

    !! candidates
    c1 = 0; c1(1) = 1
    c2 = 0; c2(1) = 1

    call ira_unify( nat1, typ1, coords1, c1, &
         nat2, typ2, coords2, c2, &
         kmax, rmat, tr, perm, dh, ierr )
    if( ierr /= 0 ) then
       call ira_get_errmsg(ierr, msg)
       write(*,*) trim(msg)
    end if
    ! write(*,*) dh
    ! write(*,*) "rmat from ira_unify"
    ! write(*,'(3f8.4)') rmat(1,:)
    ! write(*,'(3f8.4)') rmat(2,:)
    ! write(*,'(3f8.4)') rmat(3,:)
    ! write(*,*) "tr from ira_unify"
    ! write(*,'(3f8.4)') tr
    ! write(*,*)


    !! apply perm
    typ2 = typ2(perm)
    coords2 = coords2(:,perm)

    !! get svd
    call svdrot_m( nat1, typ1, coords1, nat2, typ2, coords2, rmat, tr, ierr )
    if( ierr /= 0 ) then
       write(*,*) "err svd"
       return
    end if
    ! write(*,*) "rmat from svd"
    ! write(*,'(3f8.4)') rmat(1,:)
    ! write(*,'(3f8.4)') rmat(2,:)
    ! write(*,'(3f8.4)') rmat(3,:)
    ! write(*,*) "tr from svd"
    ! write(*,'(3f8.4)') tr
    ! write(*,*)

    !! apply
    do i = 1, nat2
       coords2(:,i) = matmul(rmat, coords2(:,i)) + tr
    end do

    dh = 0.0
    !! get final dh
    do i = 1, nat1
       dx = coords2(1,i) - coords1(1,i)
       dy = coords2(2,i) - coords1(2,i)
       dz = coords2(3,i) - coords1(3,i)
       dh = max( dh, dx*dx + dy*dy + dz*dz )
    end do
    dh = sqrt(dh)

    ! write(*,*) "exit full ira"
  end function dh_full_ira


end module m_itm_core
