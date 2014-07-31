
!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!> @brief A module providing field related classes.
!>
!> @detail Both a representation of a field which provides no access to the 
!> underlying data (to be used in the algorithm layer) and an accessor class
!> (to be used in the Psy layer) are provided.


module field_mod

  use constants_mod,            only : r_def
  use function_space_mod,       only : function_space_type
  use gaussian_quadrature_mod,  only : gaussian_quadrature_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  !> Algorithm layer representation of a field.
  !>
  !> Objects of this type hold all the data of the field privately.
  !> Unpacking the data is done via the proxy type accessed by the Psy layer
  !> alone.
  !>
  type, public :: field_type
    private

    !> Each field has a pointer to the function space on which it lives
    type( function_space_type ), pointer         :: vspace => null( )
    !> Each field has a pointer to the gaussian quadrature rule which will be
    !! used to integrate over its values
    type( gaussian_quadrature_type ), pointer         &
                                      :: gaussian_quadrature => null( )
    !> Allocatable array of type real which holds the values of the field
    real(kind=r_def), allocatable         :: data( : )

  contains

    !> Function to get a proxy with public pointers to the data in a
    !! field_type.
    procedure, public :: get_proxy

    !> Sends the field contents to the log
    !! @param[in] title A title added to the log before the data is written out
    !>
    procedure, public :: print_field

    !> function returns the enumerated integer for the functions_space on which
    !! the field lives
    procedure         :: which_function_space

  end type field_type

  interface field_type

    module procedure field_constructor

  end interface

  public :: which_function_space

  !> Psy layer representation of a field.
  !>
  !> This is an accessor class that allows access to the actual field information
  !> with each element accessed via a public pointer.
  !>
  type, public :: field_proxy_type 

    private

    !> Each field has a pointer to the function space on which it lives
    type( function_space_type ), pointer, public :: vspace
    !> Each field has a pointer to the gaussian quadrature rule which will be
    !! used to integrate over its values
    type( gaussian_quadrature_type ), pointer, public &
                                      :: gaussian_quadrature
    !> Allocatable array of type real which holds the values of the field
    real(kind=r_def), public, pointer         :: data( : )

  contains
  end type field_proxy_type 

contains

  !> Function to create a proxy with access to the data in the field_type.
  !>
  !> @return The proxy type with public pointers to the elements of
  !> field_type
  type(field_proxy_type ) function get_proxy(self)
    implicit none
    class(field_type), target, intent(in)  :: self

    get_proxy % vspace                 => self % vspace
    get_proxy % gaussian_quadrature    => self % gaussian_quadrature
    get_proxy %  data                  => self % data

  end function get_proxy

  !> Construct a <code>field_type</code> object.
  !>
  !> @param [in] vector_space the function space that the field lives on
  !> @param [in] gq the gaussian quadrature rule
  !> @return self the field
  !>
  function field_constructor( vector_space, gq ) result(self)

    type(function_space_type), target, intent(in) :: vector_space
    type(gaussian_quadrature_type), optional, target, intent(in) :: gq

    type(field_type), target :: self

    self%vspace => vector_space
    if ( present( gq ) ) then
      self%gaussian_quadrature => gq
    else
      self%gaussian_quadrature => null()
    end if

    ! allocate the array in memory
    allocate(self%data(self%vspace%get_undf()))

  end function field_constructor

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------

  !> Sends the field contents to the log
  !! @param[in] title A title added to the log before the data is written out
  !>
  subroutine print_field( self, title )

    use log_mod, only : log_event, log_scratch_space, LOG_LEVEL_INFO

    implicit none

    class( field_type ), target, intent( in ) :: self

    character( * ),          intent( in ) :: title

    integer                   :: cell
    integer                   :: layer
    integer                   :: df
    integer,          pointer :: map( : )

    call log_event( title, LOG_LEVEL_INFO )

    do cell=1,self%vspace%get_ncell()
     map => self%vspace%get_cell_dofmap( cell )
      do df=1,self%vspace%get_ndf()
        do layer=0,self%vspace%get_nlayers()-1
          write( log_scratch_space, '( I4, I4, I4, F8.2 )' ) &
              cell, df, layer+1, self%data( map( df ) + layer )
          call log_event( log_scratch_space, LOG_LEVEL_INFO )
        end do
      end do
    end do

  end subroutine print_field

  function which_function_space(self) result(fs)
    implicit none
    class(field_type), intent(in) :: self
    integer :: fs
    
    fs = self%vspace%which()
    return
  end function which_function_space

end module field_mod
