module m_template

  implicit none
  private
  public :: t_template, TEMPLATE_MODE_RCUT, TEMPLATE_MODE_NN

  integer, parameter :: &
       TEMPLATE_MODE_UNKNOWN = -1, &
       TEMPLATE_MODE_RCUT = 1, &
       TEMPLATE_MODE_NN   = 2

  type :: t_template
     integer :: mode
     logical :: ignore_chem
     logical :: rescale
     !! values to retrieve structure as-input
     real :: scale
     real, dimension(3) :: origin
     !! values to store
     integer :: nat
     integer, allocatable :: typ(:)
     real, allocatable :: coords(:,:)
     integer :: hash
     real :: rcut
   contains
     procedure :: print => t_template_print
     final :: t_template_destroy
  end type t_template

  interface t_template
     procedure t_template_constructor
  end interface t_template



contains


  function t_template_constructor( nat, typ, coords_in, normalize, mode, ignore_chem ) &
       result( this )
    !! call as:
    !! type( struc_template ), pointer :: tmplt
    !!
    !! tmplt => struc_template()
    !!
    use m_itm_core, only: struc_get_scale, struc_rescale, find_dmax, order_atoms
    implicit none
    type( t_template ), pointer :: this
    integer, intent(in) :: nat
    integer, dimension(nat), intent(in) :: typ
    real, dimension(3,nat), intent(in) :: coords_in
    logical, intent(in) :: normalize
    character(*), intent(in) :: mode
    logical, intent(in) :: ignore_chem

    real, dimension(3,nat) :: coords
    real :: scale
    real, dimension(3) :: gc
    integer :: i

    !! allocate new memory
    allocate( t_template :: this )

    !! make a copy for resize
    coords = coords_in

    !! shift to gc
    gc = sum( coords, 2)/nat
    do i = 1, nat
       coords(:,i) = coords(:,i) - gc
    end do
    this% origin = gc

    !! order atoms by distance from center
    call order_atoms( nat, coords )

    !! shift such that first atom at center
    gc = coords(:,1)
    do i = 1, nat
       coords(:,i) = coords(:,i) - gc
    end do
    this% origin = this% origin + gc




    this% ignore_chem = ignore_chem
    this% rescale = normalize
    select case( mode )
       !! compare by radial cutoff
    case( "rcut" )
       this% mode = TEMPLATE_MODE_RCUT
       !! get rcut as max distance from first atom
       this% rcut = find_dmax( nat, coords )
       !! add 5 percent for some flexibility
       this% rcut = this% rcut*1.05

       !! NOTE: if we use mode=rcut, then we cannot rescale
       this% rescale = .false.

       !! compare by number of nearest neighbors
    case( "nn" )
       this% mode = TEMPLATE_MODE_NN
       !! singal that we don't know the rcut
       this% rcut = -1.0

    case default
       this% mode = TEMPLATE_MODE_UNKNOWN
       this% rcut = -999.9
    end select


    !! resize if needed
    scale = 1.0
    if( this% rescale ) then
       scale = struc_get_scale( nat, coords )
       call struc_rescale( nat, coords, 1.0/scale )
    end if
    this% scale = scale

    !! set data to template
    this% nat = nat
    allocate( this% typ, source = typ )
    allocate( this% coords, source = coords )


    !! hash
    this% hash = 0
  end function t_template_constructor


  subroutine t_template_destroy( self )
    !! deallocate memory of template
    implicit none
    type( t_template ), intent(inout) :: self

    ! write(*,*) "in t_template destroy"
    if( allocated(self% typ))deallocate( self% typ)
    if( allocated( self% coords)) deallocate( self% coords)
  end subroutine t_template_destroy


  subroutine t_template_print( self )
    !! print contents of template
    implicit none
    class( t_template ), intent(in) :: self

    integer :: i

    write(*,*) "in template_print"
    if( .not. allocated( self% typ) ) then
       write(*,*) "the template appears to have unallocated data"
       return
    end if

    write(*,*) self% nat
    write(*,'(2x,a,1x,i0,";" ,1x,a,1x,f8.4,";",1x,a,1x,l,";",1x,a,1x,l)') &
         "mode:",self% mode, &
         "rcut:",self% rcut* self% scale, &
         "ignore_chem", self% ignore_chem, &
         "rescale", self% rescale
    !! print coords non-normalized
    do i = 1, self% nat
       write(*,*) self% typ(i), self% coords(:,i)!*self% scale! + self% origin
    end do

  end subroutine t_template_print


end module m_template
