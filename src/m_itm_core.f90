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


  function struc_get_dh( nat1, typ1, coords1, nat2, typ2, coords2, kmax_factor, dthr )result(dh)
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
!       dh = max( dh, sqrt(dot_product(rdum,rdum)))
       dh = dh + dot_product( rdum, rdum )
    end do
    dh = sqrt(dh)/nat1


#ifdef DEBUG
    write(*,*) nat1+nat2
    write(*,*) 'matched'
    do i = 1, nat1
       write(*,*) 1, coords1(:,i)
       write(*,*) 2, coords2_w(:,i)
    end do
#endif




  end function struc_get_dh


  subroutine order_atoms( nat, coords )
    !! order atoms by distance from origin
    implicit none
    integer, intent(in) :: nat
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


end module m_itm_core
