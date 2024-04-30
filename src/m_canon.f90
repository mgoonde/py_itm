module m_canon

  implicit none
  private
  public :: t_canon

  !! canonical template
  type :: t_canon
     !! structure
     integer :: nat
     integer, allocatable :: typ(:)
     real, allocatable :: coords(:,:)
     integer :: pg
     !! indices
     integer :: n_similar
     integer, allocatable :: similar_canon(:)
     real, allocatable :: similar_dh(:)
   contains
     procedure :: print => t_canon_print
     final :: t_canon_destroy
  end type t_canon

  interface t_canon
     procedure t_canon_constructor
  end interface t_canon

contains

  function t_canon_constructor( nat, typ, coords ) result( this )
    use m_itm_core
    implicit none
    type( t_canon ), pointer :: this
    integer, intent(in) :: nat
    integer, intent(in) :: typ(nat)
    real,    intent(in) :: coords(3,nat)

    character(len=10) :: pg

    allocate( t_canon :: this )

    this% nat = nat
    allocate( this% typ, source=typ )
    allocate( this% coords, source=coords)

    !! get pg
    call sofi_struc_pg( nat, typ, coords, 0.3, pg, .false. )
    this% pg = pg_char2int( pg )

    !!
    this% n_similar = 0

  end function t_canon_constructor


  subroutine t_canon_destroy( this )
    implicit none
    type( t_canon ), intent(inout) :: this

    if( allocated( this% typ ))deallocate( this% typ )
    if( allocated( this% coords ))deallocate( this% coords )
    if( allocated( this% similar_canon))deallocate( this% similar_canon )
    if( allocated( this% similar_dh))deallocate( this% similar_dh )
  end subroutine t_canon_destroy


  subroutine t_canon_print( this )
    use m_itm_core, only: pg_int2char
    implicit none
    class( t_canon ), intent(in) :: this

    integer :: i
    character(len=5) :: pg

    write(*,*) "in t_canon_print"
    if( .not. allocated( this% typ )) then
       write(*,*) "t_canon appears to have unallocated data!"
       return
    end if

    pg = pg_int2char( this% pg )
    write(*,*) this% nat
    write(*,*) "PG:", trim(pg)
    do i = 1, this% nat
       write(*,*) this% typ(i), this% coords(:,i)
    end do

    write(*,*) "number of similar canons:", this% n_similar
    do i = 1, this% n_similar
       write(*,*) this% similar_canon(i), this% similar_dh(i)
    end do

  end subroutine t_canon_print

end module m_canon
