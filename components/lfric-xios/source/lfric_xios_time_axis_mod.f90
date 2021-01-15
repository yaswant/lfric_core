!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief A module containing a time axis object
!>
!> @details Fields need to be updated at different times and frequencies. The
!>          time axis object can be linked to a field to provide information
!>          of how the field should be updated with time.
module lfric_xios_time_axis_mod

  use constants_mod,        only: i_def, str_def, r_def, l_def, dp_xios
  use field_mod,            only: field_type, field_proxy_type
  use field_parent_mod,     only: field_parent_type
  use field_collection_mod, only: field_collection_type, field_collection_iterator_type
  use fs_continuity_mod,    only: W3
  use log_mod,              only: log_event, log_scratch_space, LOG_LEVEL_ERROR, LOG_LEVEL_INFO
  use linked_list_data_mod, only: linked_list_data_type

  implicit none

  private

  !> Time axis object type
  type, extends(linked_list_data_type), public :: time_axis_type

    private

    !> Name of the time_axis.
    character(str_def) :: name = 'unset'
    !> The data values of the time axis
    real(kind=r_def), allocatable :: time_data( : )
    !> The indices of the data points of the time axis
    integer(kind=i_def), allocatable :: index_data( : )
    !> The collection of fields associated with the time axis
    type(field_collection_type) :: fields
    !> Flag determining the cyclical nature of time axis
    logical(l_def) :: cyclic = .true.
    !> String identifier for units of the time axis
    character(str_def) :: time_units = 'seconds'

    procedure(update_interface), nopass, pointer :: update_behaviour => null()

  contains
    !> Initialiser for the time axis object
    procedure, public :: initialise
    !> Set the routine that updates the time axis field data
    procedure, public :: set_update_behaviour
    !> Getter for time axis name
    procedure, public :: get_name
    !> Getter for time units
    procedure, public :: get_units
    !> Procedure for cycling through time data to find the correct entry
    procedure, public :: shift_forward
    !> Returns start/end times for active time window
    procedure, public :: get_time_window
    !> Returns start/end indices for active time window
    procedure, public :: get_index_window
    !> Aligns the active time window with the current model time
    procedure, public :: align
    !> Adds a field to the associated time axis object
    procedure, public :: add_field
    !> Reads the fields associated with the time axis
    procedure, public :: update_fields
    !> Populates the model fields linked with the time axis data
    procedure, public :: populate_model_fields
  end type time_axis_type

  abstract interface

  !> Interface for updating the field data associated with the time axis
  subroutine update_interface(field_name, field_proxy, times)
    import i_def, field_proxy_type
    character(len=*),        intent(in) :: field_name
    type(field_proxy_type),  intent(inout) :: field_proxy
    integer(i_def),          intent(in) :: times(:)
  end subroutine update_interface

end interface

public :: update_interface

contains

  !> Initialise a <code>time_axis_type</code> object.
  !>
  !> @param [in] input_data The input time data array
  !> @param [in] input_index_data The corresponding indices for the data
  !> @param [in] name The time axis name
  !> @param [in] input_units The units of time data
  !> @param [in] cyclic Flag determining if axis is cyclic
  subroutine initialise(self, input_data, input_index_data, name, &
                        input_units, cyclic)

    implicit none

    class(time_axis_type),        intent(inout) :: self
    real(r_def),                  intent(in)    :: input_data(:)
    integer(i_def),               intent(in)    :: input_index_data(:)
    character(*),                 intent(in)    :: name
    character(*),   optional,     intent(in)    :: input_units
    logical(l_def), optional,     intent(in)    :: cyclic

    type(field_collection_type) :: input_fields

    self%time_data = input_data

    self%index_data = input_index_data

    self%name = name

    input_fields = field_collection_type(name=trim(name)//'_fields')
    self%fields = input_fields

    if ( present(input_units) ) self%time_units = input_units

    if ( present(cyclic) ) self%cyclic = cyclic

  end subroutine initialise

  !> Sets routine that updates field data for time axis
  !> @param [in] update_behaviour Pointer to update routine
  subroutine set_update_behaviour(self, update_behaviour)

    implicit none

    class(time_axis_type), intent(inout) :: self
    procedure(update_interface), pointer, intent(in) :: update_behaviour

    self%update_behaviour => update_behaviour

  end subroutine set_update_behaviour

  !> Returns start and end time data for active time window
  !> @result output_name The time axis name
  function get_name(self) result(output_name)

    implicit none

    class(time_axis_type), intent(inout) :: self

    character(str_def) :: output_name

    output_name = self%name

  end function get_name

  !> Returns start and end time data for active time window
  !> @result output_unitsThe time axis units
  function get_units(self) result(output_units)

    implicit none

    class(time_axis_type), intent(inout) :: self

    character(str_def) :: output_units

    output_units = self%time_units

  end function get_units

  !> Performs a cshift on the data and index data arrays.
  subroutine shift_forward(self)

    implicit none

    class(time_axis_type), intent(inout) :: self

    self%time_data = cshift(self%time_data, 1)
    self%index_data = cshift(self%index_data, 1)

  end subroutine shift_forward

  !> Returns start and end time data for active time window
  !> @result time_window The active time window of the time axis
  function get_time_window(self) result(time_window)

    implicit none

    class(time_axis_type), intent(inout) :: self

    real(r_def) :: time_window(2)

    time_window = self%time_data(1:2)

  end function get_time_window

  !> Returns start and end time indices for active time window
  !> @result time_window_index The indices of the active time window
  function get_index_window(self) result(time_window_index)

    implicit none

    class(time_axis_type), intent(inout) :: self

    integer(i_def) :: time_window_index(2)

    time_window_index = self%index_data(1:2)

  end function get_index_window

  !> Takes model time and shifts forward through time axis so active time
  !> window is aligned with model time
  !> @param [in] input_time The current time of the model
  subroutine align(self, input_time)

    implicit none

    class(time_axis_type), intent(inout) :: self
    real(r_def),           intent(in)    :: input_time

    real(r_def)    :: time_window(2)
    integer(i_def) :: n

    ! Until final clock/calendar implementation, the time is currently in the
    ! form of a real number of days (future work will introduce time units)
    time_window = self%get_time_window()
    do n=1,size(self%time_data)
      if ( (time_window(1) <= input_time) .and. (input_time < time_window(2)) ) then
        return
      else
        call self%shift_forward()
        time_window = self%get_time_window()
      end if
    end do

    ! If we're still going then the correct time is out of the time axis bounds
    ! or between the end and start of a cyclic time data set
    if ( self%cyclic ) then
      time_window = self%get_time_window()
      do n=1, size(self%index_data)
        if ( (time_window(1) > time_window(2)) ) then
          return
        else
          call self%shift_forward()
          time_window = self%get_time_window()
        end if
      end do
    else
      write( log_scratch_space, '(A,A)' ) &
        "Model has run beyond the end of the input dataset with time axis: ", &
         self%name
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

  end subroutine align

  !> Adds a field to the time axis field collection
  !> @param[in] field Field to be added to the time axis
  subroutine add_field(self, field)

    implicit none

    class(time_axis_type), intent(inout) :: self
    type(field_type), intent(in) :: field

    call self%fields%add_field(field)

  end subroutine add_field

  !> Updates associated data fields using update routine
  subroutine update_fields(self)

    implicit none

    class(time_axis_type), intent(inout) :: self

    type( field_collection_iterator_type) :: read_iter
    class( field_parent_type ), pointer   :: field => null()

    type(field_proxy_type) :: tmp_proxy


    read_iter = self%fields%get_iterator()
    do
      if ( .not.read_iter%has_next() ) exit
      field => read_iter%next()
      select type(field)
        type is (field_type)
          tmp_proxy = field%get_proxy()
          call log_event( &
            'Reading '//trim(adjustl(field%get_name())), &
            LOG_LEVEL_INFO)
          call self%update_behaviour(field%get_name(), tmp_proxy, self%index_data(1:2))
      end select
    end do

    nullify(field)

  end subroutine update_fields

  !> Populates model fields using time axis data fields
  !> @param[in] model_fields Field collection to populate from data fields
  subroutine populate_model_fields(self, model_fields)

    implicit none

    class(time_axis_type),       intent(inout) :: self
    type(field_collection_type), intent(in)    :: model_fields

    class(field_parent_type), pointer :: field => null()
    type(field_type),         pointer :: model_field => null()
    type(field_type),         pointer :: data_field => null()

    real(r_def),        allocatable :: dat_interp(:)
    character(str_def), allocatable :: field_name
    type(field_collection_iterator_type) :: pop_iter
    type(field_proxy_type) :: data_proxy, model_proxy

    pop_iter = self%fields%get_iterator()

    do
      if ( .not. pop_iter%has_next() ) exit
      field => pop_iter%next()
      field_name = field%get_name()

      if ( model_fields%field_exists(trim(field_name)) ) then
        model_field => model_fields%get_field(trim(field_name))
        model_proxy = model_field%get_proxy()
        data_field => self%fields%get_field(trim(field_name))
        data_proxy = data_field%get_proxy()

        ! Populate science fields from data - right now no interpolation is
        ! done so field is just copied via an intermediary
        allocate(dat_interp(size(model_proxy%data)))

        ! When multi-data fields are on trunk, this will be where the
        ! interpolation happens
        dat_interp = data_proxy%data

        ! Populate model field with interpolated data
        model_proxy%data = dat_interp

        ! Deallocate interpolation array ready for next field
        deallocate(dat_interp)

      end if
    end do

    nullify(field)

  end subroutine populate_model_fields

end module lfric_xios_time_axis_mod