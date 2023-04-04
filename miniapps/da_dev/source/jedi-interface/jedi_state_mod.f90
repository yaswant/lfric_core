!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing the JEDI State emulator class.
!>
!> @details This module holds a JEDI State emulator class that includes only
!>          the functionality required by the LFRic-JEDI model interface on the
!>          LFRic-API. This includes i) read/write to a defined set of fields,
!>          ii) ability to interoperate between LFRic fields and JEDI fields
!>          (Atlas fields here), and iii) storage of model_data instance to be
!>          used for IO and model time stepping.
module jedi_state_mod

  use, intrinsic :: iso_fortran_env, only : real64
  use atlas_field_emulator_mod,      only : atlas_field_emulator_type
  use atlas_field_interface_mod,     only : atlas_field_interface_type
  use jedi_geometry_mod,             only : jedi_geometry_type
  use driver_model_data_mod,         only : model_data_type
  use jedi_state_config_mod,         only : jedi_state_config_type
  use da_dev_field_meta_mod,         only : da_dev_field_meta_type
  use field_collection_mod,          only : field_collection_type
  use log_mod,                       only : log_event,          &
                                            log_scratch_space,  &
                                            LOG_LEVEL_ERROR
  use constants_mod,                 only : i_def, l_def

  implicit none

  private

type, public :: jedi_state_type
  private

  !> These fields emulate a set of Atlas fields are and used purely for testing
  type ( atlas_field_emulator_type ), allocatable :: fields(:)

  !> An object that stores the field meta-data associated with the fields
  type( da_dev_field_meta_type )                  :: field_meta_data

  !> Interface field linking the Atlas emulator fields and LFRic fields in the
  !> model data (to do field copies)
  type( atlas_field_interface_type ), allocatable :: fields_to_model_data(:)

  !> Interface field linking the Atlas emulator fields and LFRic fields in the
  !> io_collection (to do field copies)
  type( atlas_field_interface_type ), allocatable :: fields_to_io_collection(:)

  !> Model data that stores the model_data fields to propagate
  type( model_data_type ), public                 :: model_data

  !> Field collection to perform IO
  type( field_collection_type ), public           :: io_collection

  ! Here we have date_time as an integer. It will actually be an object or
  ! string that stores time to be read. Initially stored in the configuration
  ! file:
  !
  ! date: '2018-04-14T21:00:00Z'
  ! mpas define the formatting for this as:
  ! dateTimeString = '$Y-$M-$D_$h:$m:$s'
  integer                                         :: date_time

  !> The jedi_geometry object
  type( jedi_geometry_type ), pointer             :: geometry

contains

  !> Jedi state initialiser.
  procedure :: initialise => state_initialiser_read
  procedure :: state_initialiser

  !> Initialise the model_data
  procedure, private :: initialise_model_data

  !> Setup the atlas_field_interface_type that enables copying between Atlas
  !> field emulators and the LFRic fields in the model_data
  procedure, private :: setup_interface_to_model_data

  !> Setup the atlas_field_interface_type that enables copying between Atlas
  !> field emulators and the LFRic fields in io_collection
  procedure, private :: setup_interface_to_io_collection

  !> Copy the data in the LFRic fields stored in the io_collection to the
  !> internal Atlas field emulators
  procedure, private :: from_lfric_io_collection

  !> Read model fields from file into fields
  procedure, public :: read_file

  !> @todo Write model fields to file from the fields
  !> procedure, public :: write_file

  !> Get the curent time
  procedure, public :: valid_time

  !> Update the curent time
  procedure, public :: update_time

  !> Print field
  procedure, public :: print_field

  !> Copy the data in the LFRic fields stored in the model_data to the internal
  !> Atlas field emulators
  procedure, public :: from_model_data

  !> Finalizer
  final             :: jedi_state_destructor

end type jedi_state_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!-------------------------------------------------------------------------------
! Methods required by the JEDI (OOPS) model interface
!-------------------------------------------------------------------------------

!> @brief    Initialiser via read for jedi_state_type
!>
!> @param [in] geometry The geometry object required to construct the state
!> @param [in] config   The configuration object including the required
!>                      information to construct a state and read a file to
!>                      initialise the fields
subroutine state_initialiser_read( self, geometry, config )

  implicit none

  class( jedi_state_type ), intent(inout)        :: self
  type( jedi_geometry_type ), target, intent(in) :: geometry
  type( jedi_state_config_type ), intent(in)     :: config

  call self%state_initialiser( geometry, config )

  ! Initialise the Atlas field emulators via the model_data or the
  ! io_collection
  if (config%use_nl_model) then
    ! This calls the models initialise method that populates model_data and
    ! does a copy from model_data to the fields
    call self%initialise_model_data()
  else
    ! We are not running the non-linear model so read the file directly and
    ! do a copy from io_collection to the fields
    call self%read_file( config%date_time, config%read_file_prefix )
  end if

end subroutine state_initialiser_read


!> @brief    Initialiser for jedi_state_type
!>
!> @param [in] geometry The geometry object required to construct the state
!> @param [in] config   A configuration object including the required
!>                      information to construct a state
subroutine state_initialiser( self, geometry, config )

  use da_dev_model_init_mod, only : create_da_model_data
  use da_dev_driver_mod,     only : mesh, twod_mesh
  use fs_continuity_mod,     only : W3, Wtheta

  implicit none

  class( jedi_state_type ), intent(inout)        :: self
  type( jedi_geometry_type ), target, intent(in) :: geometry
  type(jedi_state_config_type), intent(in)       :: config

  ! Local
  integer( kind=i_def ) :: n_horizontal
  integer( kind=i_def ) :: n_levels
  integer( kind=i_def ) :: n_layers
  integer( kind=i_def ) :: ivar
  integer( kind=i_def ) :: n_variables
  integer( kind=i_def ) :: fs_id
  logical( kind=l_def ) :: twod_field

  ! Setup
  self%field_meta_data = config%field_meta_data
  self%date_time=config%date_time
  self%geometry => geometry
  n_variables = self%field_meta_data%get_n_variables()

  allocate( self%fields(n_variables) )

  n_horizontal=geometry%get_n_horizontal()
  n_layers=geometry%get_n_layers()

  do ivar=1,n_variables

    fs_id      = self%field_meta_data%get_variable_function_space(ivar)
    twod_field = self%field_meta_data%get_variable_is_2d(ivar)

    select case (fs_id)
    case (W3)
      if (twod_field) then
        n_levels = 1
      else
        n_levels = n_layers
      end if

    case (Wtheta)
      n_levels = n_layers + 1

    case default
      log_scratch_space = 'The requested LFRic function space is not supported.'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )

    end select

    call self%fields(ivar)%initialise(        &
                                n_levels,     &
                                n_horizontal, &
                                self%field_meta_data%get_variable_name(ivar) )
  end do

  ! Setup the io_collection (empty ready for read/write operations)
  call self%io_collection%initialise(name = 'io_collection', table_len=100)

  ! If running the model, create model data and link to fields .
  if (config%use_nl_model) then
    call create_da_model_data(mesh, twod_mesh, self%model_data)
    call self%setup_interface_to_model_data()
  end if

end subroutine state_initialiser

!> @brief    Method that returns the curent time
!>
function valid_time(self) result(date_time)

  implicit none

  class( jedi_state_type ), intent(inout) :: self
  integer( kind=i_def )                   :: date_time

  ! Here we have date_time as an integer. It will actually be an object or
  ! string that stores time to be read. Initially stored in the configuration
  ! file:
  !
  ! date: '2018-04-14T21:00:00Z'
  ! mpas define the formating for this as:
  ! dateTimeString = '$Y-$M-$D_$h:$m:$s'

  date_time = self%date_time

end function valid_time

!> @brief    A method to initialise the model_data and copy into the Atlas
!>           field emulators
!>
subroutine initialise_model_data(self)

  use da_dev_model_init_mod, only : initialise_da_model_data

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Read the state into the LFRic model data
  call initialise_da_model_data( self%model_data )

  ! Copy model_data fields to the Atlas field emulators
  call self%from_model_data()

end subroutine initialise_model_data


!> @brief    A method to update the internal Atlas field emulators
!>
!> @param [in] date_time   The data date_time to be read
!> @param [in] file_prefix Character array that specifies the file to read from
subroutine read_file(self, date_time, file_prefix)

  use da_dev_io_update_mod, only : update_io_field_collection
  use da_dev_driver_mod,   only : mesh, twod_mesh
  use lfric_xios_read_mod, only : read_state

  implicit none

  class( jedi_state_type ), intent(inout) :: self
  integer( kind=i_def ),  intent(in)      :: date_time
  character(len=*), intent(in)            :: file_prefix

  ! Set the clock to the desired read time
  call set_clock(self, date_time)

  ! Ensure the io_collection contains the variables defined in the list
  ! stored in the type
  call update_io_field_collection(self%io_collection, mesh, twod_mesh, &
                                  self%field_meta_data)

  ! Read the state into the io_collection
  call read_state( self%io_collection, prefix=file_prefix )

  ! Copy model_data fields to the Atlas field emulators
  call self%from_lfric_io_collection()

end subroutine read_file

!> Write fields stored in this state from file
!subroutine write_file()
! TBD ...
!end subroutine write_file

!------------------------------------------------------------------------------
! Local methods to support LFRic-JEDI implementation
!------------------------------------------------------------------------------

!> @brief    Setup fields_to_model_data variable that enables copying between
!>           Atlas field emulators and the LFRic fields in the model_data
!>
subroutine setup_interface_to_model_data(self)

  use da_dev_utils_mod,      only: get_model_field
  use field_mod,             only: field_type

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Local
  integer(i_def)            :: ivar
  type(field_type), pointer :: lfric_field_ptr
  real(real64), pointer     :: atlas_data_ptr(:,:)
  integer(i_def), pointer   :: horizontal_map_ptr(:)
  integer(i_def)            :: n_variables

  n_variables = self%field_meta_data%get_n_variables()

  ! Allocate space for the interface fields
  if ( allocated( self%fields_to_model_data ) ) then
    deallocate( self%fields_to_model_data )
  endif
  allocate( self%fields_to_model_data( n_variables ) )

  ! Link the Atlas emulator fields with lfric fields
  do ivar=1, n_variables

    ! Get the required data
    call get_model_field( self%field_meta_data%get_variable_name(ivar), &
                          self%model_data%depository,                   &
                          lfric_field_ptr )

    atlas_data_ptr => self%fields(ivar)%get_data()

    call self%geometry%get_horizontal_map(horizontal_map_ptr)

    call self%fields_to_model_data(ivar)%initialise( atlas_data_ptr,     &
                                                     horizontal_map_ptr, &
                                                     lfric_field_ptr )

  end do

end subroutine setup_interface_to_model_data

!> @brief    Copy from model_data to the Atlas field emulators
!>
subroutine from_model_data(self)

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Local
  integer(i_def) :: ivar

  !> @todo Will need some sort of transform for winds and
  !>       possibly other higher order elements.
  !>
  !>       call transform_winds(model_data)

  ! copy to the Atlas emulator fields
  do ivar=1,size(self%fields_to_model_data)
    call self%fields_to_model_data(ivar)%copy_from_lfric()
  end do

end subroutine from_model_data

!> @brief    Setup fields_to_io_collection variable that enables copying
!>           between Atlas field emulators and the LFRic fields in
!>           io_collection
!>
subroutine setup_interface_to_io_collection(self)

  use da_dev_utils_mod,      only: get_model_field
  use field_mod,             only: field_type

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Local
  integer(i_def)            :: ivar
  type(field_type), pointer :: lfric_field_ptr
  real(real64), pointer     :: atlas_data_ptr(:,:)
  integer(i_def), pointer   :: horizontal_map_ptr(:)
  integer(i_def)            :: n_variables

  n_variables=self%field_meta_data%get_n_variables()

  ! Allocate space for the interface fields
  if (allocated(self%fields_to_io_collection)) then
    deallocate(self%fields_to_io_collection)
  end if
  allocate(self%fields_to_io_collection(n_variables))

  ! Link the Atlas emulator fields with lfric fields
  do ivar=1,n_variables
    ! Get the required data
    call get_model_field(self%field_meta_data%get_variable_name(ivar), &
                         self%io_collection, lfric_field_ptr)
    atlas_data_ptr => self%fields (ivar)%get_data()
    call self%geometry%get_horizontal_map(horizontal_map_ptr)
    call self%fields_to_io_collection(ivar)%initialise(atlas_data_ptr,     &
                                                       horizontal_map_ptr, &
                                                       lfric_field_ptr)
  end do

end subroutine setup_interface_to_io_collection

!> @brief    Copy from model_data to the Atlas field emulators
!>
subroutine from_lfric_io_collection(self)

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Local
  integer(i_def) :: ivar

  ! Link internal Atlas field emulators to LFRic fields
  call self%setup_interface_to_io_collection()

  ! Copy the LFRic fields in the io_collection to the Atlas field emulators
  do ivar=1,size(self%fields_to_io_collection)
    call self%fields_to_io_collection(ivar)%copy_from_lfric()
  end do

end subroutine from_lfric_io_collection

!> @brief    Update the date_time by a single time-step
!>
!> @param [in] date_time_duration_dt  Update the date_time by the specified
!>                                    duration
subroutine update_time(self, date_time_duration_dt)

  implicit none

  class( jedi_state_type ), intent(inout) :: self
  integer( kind=i_def ), intent(in)       :: date_time_duration_dt

  ! Here we have date_time as an integer. It will actually be an object or
  ! string that stores time to be read. Initially stored in the configuration
  ! file:
  !
  ! date: '2018-04-14T21:00:00Z'
  ! mpas define the formating for this as:
  ! dateTimeString = '$Y-$M-$D_$h:$m:$s'

  self%date_time = self%date_time + date_time_duration_dt

end subroutine update_time

!> @brief    Set the LFRic clock to the time specified by the input date_time
!>
!> @param [in] date_time  The date_time to be used to update the LFRic clock
subroutine set_clock(self, date_time)

  use da_dev_driver_mod, only : model_clock

  implicit none

  class( jedi_state_type ), intent(inout) :: self
  integer( kind=i_def )                   :: date_time

  ! Local
  integer :: iclock
  logical :: clock_stopped

  ! date_time stored as an integer. But we will have:
  ! date: '2018-04-14T21:00:00Z'
  ! mpas define the formating for this as:
  ! dateTimeString = '$Y-$M-$D_$h:$m:$s'

  if (self%date_time > date_time) then
    write(log_scratch_space, '(A)') "The xios clock can not go backwards."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! Tick the clock to the required time but first check that the clock is
  ! running in the case its not being ticked (i.e. date_time==self%date_time).
  ! clock_stopped=.not.clock%is_running() didnt work when I tried it - set to
  ! false initially for now.
  clock_stopped=.false.
  do iclock=1,date_time-self%date_time
    clock_stopped=.not.model_clock%tick()
  end do

  ! Check the clock is still running
  if (clock_stopped) then
    write(log_scratch_space, '(A)') &
      "State::set_clock::The LFRic clock has stopped."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  self%date_time=date_time

end subroutine set_clock

!> @brief    The jedi_state_type finalizer
!>
subroutine jedi_state_destructor(self)

  implicit none

  type(jedi_state_type), intent(inout) :: self

  self % geometry => null()
  if ( allocated(self % fields ) ) then
    deallocate(self % fields )
  end if

end subroutine jedi_state_destructor

!! For testing

!> Print the 1st point in each field
subroutine print_field(self)

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Local
  real(real64), pointer :: atlas_data_ptr(:,:)
  integer               :: ivar

  ! Printing data
  do ivar=1,self%field_meta_data%get_n_variables()
    atlas_data_ptr => self%fields (ivar)%get_data()
    print*, "print_field. ivar = ", ivar, &
            "atlas_data_ptr(1,1) = ", atlas_data_ptr(1,1)
  end do

end subroutine print_field

end module jedi_state_mod
