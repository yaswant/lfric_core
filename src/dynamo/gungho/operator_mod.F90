!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!> @brief A module providing operator related classes.
!>
!> @detail 


module operator_mod

  use constants_mod,            only : r_def
  use function_space_mod,       only : function_space_type
  use mesh_mod,                 only : mesh_type
    
  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  !> Algorithm layer representation of an operator
  !>
  !> Objects of this type hold all the data of the operator privately.
  !> Unpacking the data is done via the proxy type accessed by the Psy layer
  !> alone.
  !>
  type, public :: operator_type
    private

    !> Each operator has pointers to the two function spaces on which it is
    !! defined as a map from one to the other
    type( function_space_type ), pointer         :: fs_from => null( )
    type( function_space_type ), pointer         :: fs_to => null( )
     !> Allocatable array of type real which holds the values of the operator
    real(kind=r_def), allocatable         :: local_stencil( :, :, : )
    !> Size of the outermost dimemsion of the local_stencil array, equal to 
    !! ncell*nlayers
    integer :: ncell_3d
  contains

    !> Function to get a proxy with public pointers to the data in a
    !! operator_type.
    procedure, public :: get_proxy

    !> function returns the enumerated integer for the functions_spaces which
    !! the local stencil maps from
    procedure, public :: which_fs_from

    !> function returns the enumerated integer for the functions_spaces which
    !! the local stencil mapsto
    procedure, public :: which_fs_to

    !> Returns a pointer to the mesh on which the function spaces, used by
    !> this operator, are built
    procedure, public :: get_mesh

 end type operator_type

 interface operator_type

    module procedure operator_constructor

 end interface operator_type

!  public :: which_operator

  !> Psy layer representation of an operatpr
  !>
  !> This is an accessor class that allows access to the actual operator information
  !> with each element accessed via a public pointer.
  !>
  type, public :: operator_proxy_type 

    private

    !> Each operator has pointers to the function spaces which it lives "between"
    type( function_space_type ), pointer, public :: fs_to
    type( function_space_type ), pointer, public :: fs_from
    !> Allocatable array of type real which holds the values of the operator
    real(kind=r_def), public, pointer         :: local_stencil( :, :, : )
    !> size of the outermost dimension
    integer,public                            :: ncell_3d

  contains
 end type operator_proxy_type

!------------------------------------------------------------------------------
! Module parameters
!------------------------------------------------------------------------------

contains
  

  !> Construct an <code>operator_type</code> object.
  !>
  !> @param [in] fs_from the function space that the operator maps from
  !> @param [in] fs_to the function space that the operator maps to
  !> @return self the operator
  !>
  function operator_constructor( fs_to,fs_from ) result(self)

    type(function_space_type), target, intent(in) :: fs_to
    type(function_space_type), target, intent(in) :: fs_from


    type(operator_type), target :: self

    self%fs_to   => fs_to
    self%fs_from => fs_from

    self%ncell_3d = fs_from%get_ncell() * fs_from%get_nlayers()
    ! allocate the array in memory
    allocate(self%local_stencil( fs_to%get_ndf(),fs_from%get_ndf(), self%ncell_3d ) )
    
  end function operator_constructor

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !--------------------------------------------------------------------------

  !> Function to create a proxy with access to the data in the operator_type.
  !>
  !> @return The proxy type with public pointers to the elements of
  !> operator_type
  type(operator_proxy_type ) function get_proxy(self)
    implicit none
    class(operator_type), target, intent(in)  :: self

    get_proxy % fs_from                 => self % fs_from
    get_proxy % fs_to                   => self % fs_to
    get_proxy % local_stencil           => self % local_stencil
    get_proxy % ncell_3d                =  self%ncell_3d

  end function get_proxy


  function which_fs_from(self) result(fs)
    implicit none
    class(operator_type), intent(in) :: self
    integer :: fs
    
    fs = self%fs_from%which()

    return
  end function which_fs_from

  function which_fs_to(self) result(fs)
    implicit none
    class(operator_type), intent(in) :: self
    integer :: fs
    
    fs = self%fs_to%which()

    return
  end function which_fs_to

  function get_mesh(self) result(mesh)

    implicit none

    class (operator_type), intent(in) :: self
    type(mesh_type), pointer :: mesh

    mesh => self%fs_from%get_mesh()

    return
  end function get_mesh

end module operator_mod
