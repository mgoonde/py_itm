module m_itm_compute

  use m_itm
contains


  !> @details
  !! This routine checks if there are templates added before the computation
  !! has started, and checks each such template against all sites.
  subroutine check_pre_templates( fptr, dthr )
    use m_itm_core
    use omp_lib
    implicit none
    type( t_itm_ptr ), pointer :: fptr
    real, intent(in) :: dthr

    integer :: itemplate, isite

    !!
    !! loop over any existing templates
    !!
    do itemplate = 1, fptr% ntemplate
       !! omp do
       do isite = 1, fptr% nat
       end do
       !! omp end do
    end do

  end subroutine check_pre_templates



  !> @details
  !! Create a first-pass list of templates which are variant to rigid transformations, i.e.
  !! the distance d(A,B) is computed with one-shot CShDA.
  subroutine create_rough_list( fptr, dthr )
    use m_itm_core
    use omp_lib
    implicit none
    type( t_itm_ptr ), pointer :: fptr
    real, intent(in) :: dthr

    real :: buffer_thr
    integer :: isite, i, ierr
    integer :: nthread
    integer :: nat_loc, nat2
    integer, allocatable :: typ_loc(:), typ2(:)
    real, allocatable :: coords_loc(:,:), coords2(:,:)
    integer, allocatable :: perm(:)
    integer :: me
    real :: dh

    !! threshold for buffer zone around template
    buffer_thr = dthr/2.0

    !$OMP PARALLEL
    nthread = OMP_GET_MAX_THREADS()
    !$OMP END PARALLEL

    !!
    !! outer loop over isite which do not have a matching template:
    !! define this local conf as new template, and then inner loop over all sites
    !!
    do isite = 1, fptr% nat
       !!
       !! skip sites which are already matched
       !!
       if( fptr% site_template(isite) > 0 ) cycle
       !!
       !! for each local site:
       !! select local conf, and set it as new template
       !!
       call get_local_conf( fptr, isite, nat_loc, typ_loc, coords_loc )

       !!
       !! add template
       !!
       ierr = itm_add_template( fptr, nat_loc, typ_loc, coords_loc, .false., "nn", .false. )

       !!
       !! set data of current site:
       !! dh=0.0, and template index is the last template added to list
       !!
       fptr% site_dh( isite ) = 0.0
       fptr% site_template( isite ) = fptr% ntemplate
       fptr% perm_site2rough( isite ) = [ (i, i=1,nat_loc) ]

       !!
       !! inner loop over all isite
       !!
       !$OMP PARALLEL PRIVATE( me )
       me = OMP_GET_THREAD_NUM()
       !$OMP DO PRIVATE(i, nat2, typ2, coords2, dh, perm, ierr ) SCHEDULE( static )
       do i = 1, fptr% nat
          !! "diagonal element" is zero.
          if( i == isite ) cycle

          !!
          !! at this point, comparisons are done only with cshda, which should have radial buffer region
          !! identical in all directions (is a 'sphere') around a point, so a point which is within
          !! distance buffer_thr=dthr/2.0 cannot actually be closer to any other template than the original
          if( fptr% site_dh(i) < buffer_thr ) cycle

          !! get local conf for site=i and cmopute dh
          call get_local_conf( fptr, i, nat2, typ2, coords2 )
          allocate( perm(1:nat2))
          dh = dh_cshda( nat_loc, typ_loc, coords_loc, nat2, typ2, coords2, dthr, perm )

          !! distance is below thr, there are two cases to consider:
          !! 1) the site i is not yet matched:
          !!      accept the current template
          !! 2) the site is already matched:
          !!      accept the current template only if this dh is lower than previous
          if( dh < dthr ) then
             select case( fptr% site_template(i) )
             case( -1 )
                !! site is not matched yet, accept template
                fptr% site_template(i) = fptr% ntemplate
                fptr% site_dh(i) = dh
                fptr% perm_site2rough(i) = perm

             case default
                !! site is already matched, accept only if dh is lower than previous match
                if( dh < fptr% site_dh(i) ) then
                   fptr% site_template(i) = fptr% ntemplate
                   fptr% site_dh(i) = dh
                   fptr% perm_site2rough(i) = perm
                end if

             end select
             !!
          end if
          !!
          deallocate( perm )
          deallocate( typ2, coords2 )

       end do
       !$OMP END DO
       !$OMP END PARALLEL
       !!
       deallocate( typ_loc, coords_loc )
    end do

  end subroutine create_rough_list



  !> @details
  !! Create a list of canonical templates: a canonical template is invariant to rigid
  !! transformations, i.e. the d(A,B) is computed by full IRA.
  subroutine create_canon_list( fptr, dthr )
    use m_itm_core
    use m_template
    use m_canon
    use omp_lib
    implicit none
    type( t_itm_ptr ), pointer :: fptr
    real, intent(in) :: dthr

    integer :: i, j, ierr
    real :: dh
    type( t_template ), pointer :: t1, t2
    type( t_canon ), pointer :: c
    integer, allocatable :: perm(:)

    write(*,*) "got ntemplate:", fptr% ntemplate

    ! fptr% ntemplate=2
    do i = 1, fptr% ntemplate
       !!
       !! get template
       !!
       nullify( t1 )
       t1 => fptr% get_template(i)


       !! this template already has canon, can skip it
       if( t1% canon_idx > 0 ) cycle
       ! if( t1% canon_dh < dthr/2.0 ) cycle

       !! create new canonical template
       ierr = itm_add_canon( fptr, t1% nat, t1% typ, t1% coords, dthr )

       t1% is_canon = 1
       t1% canon_idx = fptr% ncanon
       t1% canon_dh = 0.0
       allocate(t1% perm_rough2canon, source=[(j,j=1,t1%nat)])

       !$OMP PARALLEL PRIVATE( j, t2, perm, dh )
       !$OMP DO SCHEDULE( dynamic, 5 )
       do j = i+1, fptr% ntemplate
          !!
          if( i == j ) cycle
          !!
          !! get template
          !!
          nullify( t2 )
          t2 => fptr% get_template(j)

          !! template already linked to some canon, with dh < d_buffer
          if( t2% canon_dh < dthr/2.0 ) cycle

          !! compare the templates with full IRA
          allocate( perm(1:t2% nat))
          dh = dh_full_ira( t1% nat, t1% typ, t1% coords, t2% nat, t2% typ, t2% coords, 1.2, perm )

          ! if( i==1 .and. j==2) then
          !    write(*,*) "1,2", dh
          !    write(*,*) t2% canon_idx
          !    write(*,*) t2% canon_dh
          !    stop
          ! end if

          !! dh is less than dthr, two cases:
          !! 1) tmplt does not have a canon assignment yet, accept this
          !! 2) tmplt already has a canon:
          !!        overwrite with current only if dh < canon_dh
          if( dh < dthr ) then

             write(*,*) "compare:",i,j,dh
             select case( t2% canon_idx )
             case( -1 )
                t2% canon_idx = fptr% ncanon
                t2% canon_dh = dh
                allocate( t2% perm_rough2canon, source=perm )

             case default
                if( dh < t2% canon_dh ) then
                   t2% canon_idx = fptr% ncanon
                   t2% canon_dh = dh
                   if( allocated(t2% perm_rough2canon))deallocate( t2% perm_rough2canon )
                   allocate( t2% perm_rough2canon, source=perm )
                end if

             end select
          end if


          !! keep perm_rough2canon


          deallocate( perm )
       end do
       !$OMP END DO
       !$OMP END PARALLEL
    end do

    !!
    !!=======
    !! NOTE: Need to "symmetrize" the similar_canon
    !!=======
    !!




    ! do i = 1, fptr% ntemplate
    !    nullify(t1)
    !    t1 => fptr% get_template(i)
    !    write(*,*) i, t1% canon_idx
    ! end do

  end subroutine create_canon_list



  !> @details
  !! Assign the canon templates to all sites
  subroutine assign_canon2site( fptr, dthr )
    use m_template
    use m_itm_core
    implicit none
    type( t_itm_ptr ), pointer :: fptr
    real, intent(in) :: dthr

    integer :: i, tmplt_idx, canon_idx, j
    type( t_template ), pointer :: t
    type( t_canon ), pointer :: c, c_sim

    integer :: ncanon

    integer :: nat_loc
    integer, allocatable :: typ_loc(:)
    real, allocatable :: coords_loc(:,:)
    integer, allocatable :: perm(:)

    real :: dd, dh
    real :: rmat(3,3), tr(3)
    integer :: ierr

    !! omp do
    !$OMP PARALLEL PRIVATE( i, nat_loc, typ_loc, coords_loc, tmplt_idx, t, c, j, &
    !$OMP                   canon_idx, dd, perm, rmat, tr, dh, c_sim )
    !$OMP DO SCHEDULE( static, 500 )
    do i = 1, fptr% nat
    ! do i = 1, 2

       !! get local conf
       call get_local_conf( fptr, i, nat_loc, typ_loc, coords_loc )

       !! get tmplt
       tmplt_idx = fptr% site_template(i)
       nullify(t)
       t => fptr% get_template(tmplt_idx)

       !! permute local2template
       ! allocate( perm, source=fptr% perm_site2rough(i)% perm)
       perm = fptr% perm_site2rough(i)
       typ_loc = typ_loc(perm)
       coords_loc = coords_loc(:,perm)

       !! permute to template local order
       typ_loc = typ_loc(t%order)
       coords_loc = coords_loc(:,t%order)

       ! write(*,*) nat_loc
       ! write(*,*) "loc", i
       ! do j = 1, nat_loc
       !    write(*,*) typ_loc(j), coords_loc(:,j)
       ! end do
       ! write(*,*) nat2
       ! write(*,*) "tmplt", tmplt_idx
       ! do j = 1, nat2
       !    ! write(*,*) t% typ(j), t% coords(:,j)
       !    write(*,*) typ2( j ), coords2(:,j)
       ! end do

       dd = 0.0
       do j = 1, nat_loc
          ! write(*,*) i, norm2( coords_loc(:,j) - coords2(:,j) )
          dd = max( dd, norm2(coords_loc(:,j) - t% coords(:,j) ))
       end do
       ! if( dd.gt.1e-1)write(*,*) "loc tmpl",i, tmplt_idx, dd


       !! permute template2canon
       typ_loc = typ_loc( t% perm_rough2canon )
       coords_loc = coords_loc(:, t% perm_rough2canon )

       !! get the canon
       canon_idx = t% canon_idx
       nullify(c)
       c => fptr% get_canon( canon_idx )


       !! first guess is current canon:
       !! one-shot svd+cshda for local and canon
       call svdrot_m( c% nat, c% typ, c% coords, nat_loc, typ_loc, coords_loc, rmat, tr, ierr )
       do j = 1, nat_loc
          coords_loc(:,j) = matmul(rmat, coords_loc(:,j)) + tr
       end do
       dd = dh_cshda( c%nat, c%typ, c%coords, nat_loc, typ_loc, coords_loc, 999.9, perm )
       ! if(dd.gt.1e-2)write(*,*) i, canon_idx, dd


       !! worst case... this should be check on d > (dAB - dthr )
       if( dd > dthr/2.0 ) then
          !! loop through all similar to canon
          !! to see if any has lower dh
          do j = 1, c% n_similar
             !! full IRA dh
             nullify( c_sim )
             c_sim => fptr% get_canon( c% similar_canon(j) )
             !! get dh
             dh = dh_full_ira( c_sim%nat, c_sim%typ, c_sim%coords, nat_loc, typ_loc, coords_loc, 1.4, perm )
             if( dh < dd ) then
                write(*,*) "lower dh:",i,canon_idx,c%similar_canon(j),dd, dh
                canon_idx = c%similar_canon(j)
                dd = dh
             end if
          end do
       endif

       deallocate( perm )
       ! deallocate( typ_loc, coords_loc )

       !! set final data
       fptr% site_template(i) = canon_idx
       fptr% site_dh(i) = dd

       !! template is canon, do nothing
       ! if( t% is_canon > 0 ) cycle

       !! template is not canon, replace
       ! write(*,*) i, tmplt_idx, t% canon_idx, fptr% site_dh(i), t% canon_dh
       ! fptr% site_template(i) = map_canon2idx( t% canon_idx )
       ! fptr% site_template(i) = t% canon_idx
       ! fptr% site_dh(i) = fptr% site_dh(i) + t% canon_dh
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    !! omp end do




  end subroutine assign_canon2site



end module m_itm_compute

