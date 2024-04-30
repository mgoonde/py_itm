module m_perm

  !! define the permutation type and functions to deal with it
  type :: t_perm
     integer, allocatable :: perm(:)
   contains
     final :: t_perm_destroy
  end type t_perm


  interface assignment(=)
     module procedure assign_perm_tperm, assign_tperm_perm
  end interface assignment(=)


contains


  !! destructor
  subroutine t_perm_destroy( self )
    type( t_perm ), intent(inout) :: self
    if( allocated( self% perm))deallocate( self% perm )
  end subroutine t_perm_destroy

  !!
  !! assignment for obtaining perm array from t_perm:
  !!
  !! integer, allocateble :: lhs(:)
  !! type( t_perm ) :: rhs
  !!
  !! lhs = rhs
  !!
  subroutine assign_perm_tperm( lhs, rhs )
    integer, allocatable, intent(out) :: lhs(:)
    class( t_perm ), intent(in) :: rhs
    allocate( lhs, source=rhs% perm )
  end subroutine assign_perm_tperm
  !!
  !! assignment for setting t_perm%perm(:) array from perm(:)
  !!
  !! type( t_perm ) :: lhs
  !! integer, dimension(n) :: rhs
  !!
  !! lhs = rhs
  !!
  subroutine assign_tperm_perm( lhs, rhs )
    class( t_perm ), intent(inout) :: lhs
    integer, intent(in) :: rhs(:)
    !! test if t_perm already filled, if yes overwrite
    if( allocated( lhs% perm ))deallocate(lhs% perm)
    allocate( lhs% perm, source = rhs)
  end subroutine assign_tperm_perm


end module m_perm
