!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a field meta-data type.
!>
!> @details Defines an object to store the meta-data associated with
!>          collections of fields used by the JEDI model-interface.
!>
module da_dev_field_meta_mod

  use constants_mod,   only : i_def, str_def, l_def, i_native
  use log_mod,         only : log_event, log_scratch_space, LOG_LEVEL_ERROR

  implicit none

  !private

  type, public :: da_dev_field_meta_type

    !private

    !> The number of variables
    integer( kind=i_def )                 :: n_variables

    !> A list of variable field names
    character( len=str_def ), allocatable :: variable_names(:)

    !> A list of varible field function_spaces
    integer( kind=i_def ), allocatable    :: variable_function_spaces(:)

    !> A list of logicals that defines if the fields are 2D
    logical( kind=l_def ), allocatable    :: variable_is_2d(:)

  contains

    !> Initialiser
    procedure :: initialise

    !> Copy constructor
    procedure, public :: da_dev_field_meta_type_assign

    !> Get the number of variables
    procedure :: get_n_variables

    !> Get the variable name for a given index
    procedure :: get_variable_name

    !> Get the variable functionspace for a given index
    procedure :: get_variable_function_space

    !> Get the variable logical defining if field is 2d for a given index
    procedure :: get_variable_is_2d

    !> Check the index is within bounds
    procedure :: check_index

    !> da_dev_field_meta_type finalizer
    final     :: da_dev_field_meta_destructor

    !> Override default assignment operator.
    generic :: assignment(=) => da_dev_field_meta_type_assign

  end type da_dev_field_meta_type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains


!> @brief    Initialiser for da_dev_field_meta_type
!>
!> @param[in] variable_function_spaces a list of function_spaces
!> @param[in] variable_is_2d           a list of logicals defining if the field
!>                                     uses the specified 2D mesh
!> @param[in] variable_names           a list of required variable names
subroutine initialise( self, variable_names, variable_function_spaces, &
                       variable_is_2d )

  implicit none

  class( da_dev_field_meta_type ), intent(inout) :: self
  character(*), intent(in)                       :: variable_names(:)
  integer(kind=i_native), intent(in)             :: variable_function_spaces(:)
  logical(kind=l_def), intent(in)                :: variable_is_2d(:)

  ! Get and check the size of the input data is consistent
  self%n_variables=size(variable_names)
  if ( self%n_variables /= size(variable_function_spaces ) .or. &
       self%n_variables /= size(variable_is_2d) ) then
    write(log_scratch_space, '(A)') &
      "Check the size of the input arrays, they should be consistent."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! Allocate and copy input meta-data
  allocate( self%variable_names(self%n_variables) )
  allocate( self%variable_function_spaces(self%n_variables) )
  allocate( self%variable_is_2d(self%n_variables) )

  self%variable_names = variable_names
  self%variable_function_spaces = variable_function_spaces
  self%variable_is_2d = variable_is_2d

end subroutine initialise

!> @brief    Get the number of variables
!>
!> @return n_variables The number of variables
function get_n_variables(self) result(n_variables)

  implicit none

  class( da_dev_field_meta_type ) :: self
  integer( kind=i_def )           :: n_variables
  n_variables = self%n_variables

end function get_n_variables

!> @brief    Get the variable name for a given index
!>
!> @param[in] variable_index the index in the array to return
!> @return    variable_name  the variable name stored in requested index
function get_variable_name(self, variable_index) result(variable_name)

  implicit none

  class( da_dev_field_meta_type ), target :: self
  integer( kind=i_def )                   :: variable_index
  character( len=str_def ), pointer       :: variable_name

  call self%check_index(variable_index,"get_variable_name")
  variable_name => self%variable_names(variable_index)

end function get_variable_name

!> @brief    Get the variable functionspace for a given index
!>
!> @param[in] variable_index          the index in the array to return
!> @return    variable_function_space the variable function space stored in
!>                                    requested index
function get_variable_function_space(self, variable_index) &
           result(variable_function_space)

  implicit none

  class( da_dev_field_meta_type ) :: self
  integer( kind=i_def )           :: variable_index
  integer( kind=i_def )           :: variable_function_space

  call self%check_index(variable_index,"get_variable_function_space")
  variable_function_space = self%variable_function_spaces(variable_index)

end function get_variable_function_space

!> @brief    Get the variable logical defining if field is 2D for a given index
!>
!> @param[in] variable_index the index in the array to return
!> @return    variable_is_2d the logical defining if field is 2D for the
!>                           requested index
function get_variable_is_2d(self, variable_index) result(variable_is_2d)

  implicit none

  class( da_dev_field_meta_type ) :: self
  integer( kind=i_def )           :: variable_index
  logical( kind=l_def )           :: variable_is_2d

  call self%check_index(variable_index,"get_variable_is_2d")
  variable_is_2d = self%variable_is_2d(variable_index)

end function get_variable_is_2d

!> @brief    Check the index is within bounds
!>
!> @param[in] variable_index the index to check
!> @param[in] calling_routine a character to print if test fails
subroutine check_index(self, variable_index, calling_routine)

  implicit none

  class( da_dev_field_meta_type ), intent(in) :: self
  integer( kind=i_def ), intent(in)           :: variable_index
  character( len=* ), intent(in)              :: calling_routine

  if (variable_index < 1 .or. variable_index > self%n_variables) then
    write(log_scratch_space, '(2A,I0,A)') trim(calling_routine), &
      ": The variable_index (=", variable_index, ") which is out of bounds."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

end subroutine check_index

!> @brief Copy constructor.
!>
!> @param[out] dest Object to receive copy.
!> @param[in] source Object to be copied.
!>
subroutine da_dev_field_meta_type_assign(dest, source)

  implicit none

  class( da_dev_field_meta_type ), intent(out) :: dest
  class( da_dev_field_meta_type ), intent(in)  :: source


  ! Copy the size
  dest%n_variables=source%n_variables

  ! Deallocate if already allocated
  if ( allocated(dest%variable_names) ) deallocate(dest%variable_names)
  if ( allocated(dest%variable_is_2d) ) deallocate(dest%variable_is_2d)
  if ( allocated(dest%variable_function_spaces) ) &
    deallocate(dest%variable_function_spaces)

  ! Allocate
  allocate( dest%variable_names(dest%n_variables) )
  allocate( dest%variable_is_2d(dest%n_variables) )
  allocate( dest%variable_function_spaces(dest%n_variables) )

  ! Copy input meta-data
  dest%variable_names = source%variable_names
  dest%variable_is_2d = source%variable_is_2d
  dest%variable_function_spaces = source%variable_function_spaces

end subroutine da_dev_field_meta_type_assign

!> @brief    da_dev_field_meta_type finalizer
!>
subroutine da_dev_field_meta_destructor(self)

  implicit none

  type( da_dev_field_meta_type ), intent(inout)    :: self

  if ( allocated(self%variable_names) ) deallocate(self%variable_names)
  if ( allocated(self%variable_function_spaces) ) &
                                      deallocate(self%variable_function_spaces)
  if ( allocated(self%variable_is_2d) ) deallocate(self%variable_is_2d)

end subroutine da_dev_field_meta_destructor

end module da_dev_field_meta_mod
