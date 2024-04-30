module m_itm

  !! this module defines the t_itm_ptr, and its internal routines


  use m_linked_list
  use m_itm_core
  use m_template
  use m_perm
  implicit none

  real, parameter :: kmax_factor = 1.8

  type :: t_itm_ptr
     integer :: ntemplate
     type( linked_list ) :: template_list

     integer :: nfast
     type( linked_list ) :: fast_list

     real :: dthr
     integer :: nat                         !! number of all atoms
     integer, allocatable :: neighlist(:,:) !! neighbor list from ovito
     real, allocatable :: veclist(:,:)      !! vector list from ovito
     integer, allocatable :: count_n(:)     !! cumulative sum of neighs
     integer, allocatable :: typ(:)         !! all atomic types
     !!
     !! results
     integer, allocatable :: site_template(:)
     real, allocatable :: site_dh(:)
     integer, allocatable :: site_pg(:)
     real, allocatable :: site_strain(:)

     type( t_perm ), allocatable :: perm_site2rough(:)
   contains
     procedure :: get_template
  end type t_itm_ptr



contains


  function itm_init()result(fptr)
    implicit none
    type( t_itm_ptr ), pointer :: fptr

    allocate( t_itm_ptr :: fptr )

    !! initialize linked list
    fptr% template_list = linked_list()
    fptr% ntemplate = 0

    fptr% nfast = 0
    fptr% fast_list = linked_list()

    fptr% nat = -1

    !! initialize to large
    fptr% dthr = 999.9

  end function itm_init

  subroutine itm_close(fptr)
    implicit none
    type( t_itm_ptr ), pointer :: fptr

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

    if( allocated( fptr% perm_site2rough)) deallocate( fptr% perm_site2rough)
    deallocate( fptr )

  end subroutine itm_close


  function itm_add_template( fptr, nat, typ, coords, normalize, mode, ignore_chem )result(ierr)
    use m_itm_core, only: dh_cshda
    implicit none
    type( t_itm_ptr ), pointer :: fptr
    integer, intent(in) :: nat
    integer, intent(in) :: typ(nat)
    real, intent(in) :: coords(3,nat)
    logical, intent(in) :: normalize
    character(*), intent(in) :: mode
    logical, intent(in) :: ignore_chem
    integer :: ierr

    type( t_template ), pointer :: tmplt, tmplt_1
    class(*), pointer :: p
    integer, allocatable :: perm(:)
    integer :: i
    real :: dh

    !! create new template
    tmplt => t_template( nat, typ, coords, normalize, mode, ignore_chem )

    !! create comparison to all previous templates in list
    !!
    !! define all previous templates in list as similar (compute dh)
    allocate( tmplt% similar_template_idx(1:fptr% ntemplate) )
    allocate( tmplt% similar_template_dh(1:fptr% ntemplate) )
    do i = 1, fptr% ntemplate

       tmplt% similar_template_idx(i) = i

       !! get template i
       nullify( tmplt_1 )
       tmplt_1 => get_template( fptr, i )
       !! compare with cshda
       allocate( perm(1:tmplt_1% nat) )
       dh = dh_cshda( tmplt% nat, tmplt% typ, tmplt% coords, &
                      tmplt_1% nat, tmplt_1% typ, tmplt_1% coords, fptr% dthr, perm )
       deallocate( perm )

       tmplt% similar_template_dh(i) = dh
    end do

    !! add template pointer to list
    fptr% ntemplate = fptr% ntemplate + 1
    nullify(p)
    p => tmplt
    call fptr% template_list% add_pointer( fptr% ntemplate, p )


    ierr = 0

  end function itm_add_template



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



  subroutine get_local_conf( fptr, isite, nat_loc, typ_loc, coords_loc )
    !! helper routine to obtain local conf around site isite
    use m_itm_core
    implicit none
    type( t_itm_ptr ), pointer :: fptr
    integer, intent(in) :: isite
    integer, intent(out) :: nat_loc
    integer, allocatable, intent(out) :: typ_loc(:)
    real, allocatable, intent(out) :: coords_loc(:,:)

    integer, allocatable :: idx(:)
    integer :: i

    call extract_elements( &
         size(fptr% neighlist, 1), size( fptr% neighlist, 2), fptr% neighlist,&
         size(fptr% veclist, 2), fptr% veclist, fptr% nat, fptr% count_n,&
         isite, nat_loc, idx, coords_loc )

    allocate( typ_loc(1:nat_loc) )
    do i = 1, nat_loc
       typ_loc(i) = fptr% typ( idx(i))
    end do

    deallocate( idx )
  end subroutine get_local_conf



end module m_itm
