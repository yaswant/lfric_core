!-------------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief XIOS wrapper for use in diagnostics
!> @details XIOS wrapper to be used in the initialisation of diagnostic fields

module lfric_xios_diag_mod

  use constants_mod,                  only: l_def, i_def, str_def
#ifdef UNIT_TEST
  use lfric_xios_mock_mod,            only:                                   &
                                        xios_is_valid_file,                   &
                                        xios_is_defined_file_attr,            &
                                        xios_get_file_attr,                   &
                                        xios_is_valid_field,                  &
                                        xios_is_defined_field_attr,           &
                                        xios_get_field_attr,                  &
                                        xios_field_is_active,                 &
                                        xios_get_axis_attr,                   &
                                        lfric_xios_mock_pull_in
#else
  use lfric_xios_mock_mod,            only: lfric_xios_mock_pull_in
  use xios,                           only:                                   &
                                        xios_is_valid_file,                   &
                                        xios_is_defined_file_attr,            &
                                        xios_get_file_attr,                   &
                                        xios_is_valid_field,                  &
                                        xios_is_defined_field_attr,           &
                                        xios_get_field_attr,                  &
                                        xios_field_is_active,                 &
                                        xios_get_axis_attr

#endif

  use log_mod,                         only:                                  &
                                        log_event,                            &
                                        log_level_error,                      &
                                        log_level_warning,                    &
                                        log_scratch_space

  implicit none

  private

  public ::                                                                    &
    file_is_enabled,                                                           &
    file_is_tagged,                                                            &
    get_file_name,                                                             &
    field_is_valid,                                                            &
    field_is_enabled,                                                          &
    field_is_active,                                                           &
    get_field_order,                                                           &
    get_field_grid_ref,                                                        &
    get_field_domain_ref,                                                      &
    get_field_axis_ref,                                                        &
    get_axis_dimension

contains

  !> @brief Return true if an only if an XIOS output file is enabled.
  !> @param[in]    unique_id  XIOS id of the file
  !> @return                  Boolean representing enabled/disabled status
  function file_is_enabled(file_id) result(enabled)
    implicit none
    character(*), intent(in) :: file_id
    logical(l_def) :: has_enabled_flag
    logical(l_def) :: enabled
    if (.not. xios_is_valid_file(file_id)) then
      write(log_scratch_space, '(A, A)')                                      &
        'Invalid XIOS file:', file_id
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if
    call xios_is_defined_file_attr(file_id, enabled=has_enabled_flag)
    if (has_enabled_flag) then
       call xios_get_file_attr(file_id, enabled=enabled)
       ! explicitly enabled/disabled
    else
       enabled = .true. ! enabled by default
    end if
  end function file_is_enabled

  !> @brief Scan comment field of a file for (boolean) tag, e.g, [[diag]]
  !> @param[in]    file_id      XIOS id of the field
  !> @param[in]    tag          Tag string, e.g., "diag"
  !> @return                    True if and only if the file is tagged
  function file_is_tagged(file_id, tag) result(tagged)
    implicit none
    character(len=*), intent(in) :: file_id
    character(len=*), intent(in) :: tag
    logical(l_def) :: tagged
    logical(l_def) :: def
    character(str_def) :: comment
    character(str_def) :: pattern
    if (.not. xios_is_valid_file(file_id)) then
      write(log_scratch_space, '(A, A)')                                      &
        'Invalid XIOS file:', file_id
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if
    tagged = .false.
    call xios_is_defined_file_attr(file_id, comment=def)
    if (def) then
      call xios_get_file_attr(file_id, comment=comment)
      pattern = '[[' // tag //']]'
      tagged = (index(comment, pattern) /= 0)
    end if
  end function file_is_tagged

  !> @brief Return the name of an XIOS file
  !> @param[in]    file_id      XIOS id of the field
  !> @return                    The file name
  function get_file_name(file_id) result(name)
    implicit none
    character(len=*), intent(in) :: file_id
    logical(l_def) :: def
    character(str_def) :: name
    if (.not. xios_is_valid_file(file_id)) then
      write(log_scratch_space, '(A, A)')                                      &
        'Invalid XIOS file:', file_id
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if
    call xios_is_defined_file_attr(file_id, name=def)
    if (def) then
      call xios_get_file_attr(file_id, name=name)
    else
      name = file_id
    end if
  end function get_file_name

  !> @brief Return true if an only if an XIOS field id is valid.
  !> @param[in]    unique_id    XIOS id of the field
  !> @return                    Boolean representing validity
  function field_is_valid(unique_id) result(valid)
    implicit none
    character(*), intent(in) :: unique_id
    logical(l_def) :: valid
    valid = xios_is_valid_field(unique_id)
  end function field_is_valid

  !> @brief Return true if an only if an XIOS output field is enabled.
  !> @param[in]              unique_id    XIOS id of the field
  !> @return                    Boolean representing enabled/disabled status
  function field_is_enabled(unique_id) result(enabled)
    implicit none
    character(*), intent(in) :: unique_id
    logical(l_def) :: has_enabled_flag
    logical(l_def) :: enabled
    if (.not. xios_is_valid_field(unique_id)) then
      write(log_scratch_space, '(A, A)')                                      &
      'Invalid XIOS field:', unique_id
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    else
      call xios_is_defined_field_attr(unique_id, enabled=has_enabled_flag)
      if (has_enabled_flag) then
        call xios_get_field_attr(unique_id, enabled=enabled)
        ! explicitly enabled/disabled
      else
        enabled = .true. ! enabled by default
      end if
    end if
  end function field_is_enabled

  !> @brief Return true if and only if an XIOS field is active.
  !> @param[in]              unique_id            XIOS id of the field
  !> @param[in]              at_current_timestep  Active at current timestep?
  !> @return                 Boolean representing active/inactive status
  function field_is_active(unique_id, at_current_timestep) result(active)
    implicit none
    character(*), intent(in) :: unique_id
    logical(l_def) :: at_current_timestep
    logical(l_def) :: active
    if (.not. xios_is_valid_field(unique_id)) then
      write(log_scratch_space, '(A, A)')                                      &
      'Invalid XIOS field:', unique_id
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    else
      active = xios_field_is_active(unique_id, at_current_timestep)
    end if
  end function field_is_active

  !> @brief Return the interpolation order of a field.
  !> @param[in]    unique_id    XIOS id of the field
  !> @return                    The interpolation order
  function get_field_order(unique_id) result(order)
   implicit none
   character(*), intent(in) :: unique_id
   integer(i_def) :: order
   order = 0 ! FOR NOW, to be extracted from the field comment?
  end function get_field_order

  !> @brief Return the XIOS grid reference of a field.
  !> @param[in]    unique_id    XIOS id of the field
  !> @return                    The grid reference
  function get_field_grid_ref(unique_id) result(grid_ref)
    implicit none
    character(*), intent(in) :: unique_id
    character(str_def) :: grid_ref
    logical(l_def) :: has_grid
    call xios_is_defined_field_attr(unique_id, grid_ref=has_grid)
    if (has_grid) then
       call xios_get_field_attr(unique_id, grid_ref=grid_ref)
    else
       grid_ref = ""
    end if
  end function get_field_grid_ref

  !> @brief Return the XIOS domain reference of a field.
  !> @param[in]    unique_id    XIOS id of the field
  !> @return                    The domain reference
  function get_field_domain_ref(unique_id) result(domain_ref)
    implicit none
    character(*), intent(in) :: unique_id
    character(str_def) :: domain_ref
    logical(l_def) :: has_domain
    call xios_is_defined_field_attr(unique_id, domain_ref=has_domain)
    if (has_domain) then
       call xios_get_field_attr(unique_id, domain_ref=domain_ref)
    else
       domain_ref = ""
    end if
  end function get_field_domain_ref

  !> @brief Return the XIOS axis reference of a field.
  !> @param[in]    unique_id    XIOS id of the field
  !> @return                    The axis reference
  function get_field_axis_ref(unique_id) result(axis_ref)
    implicit none
    character(*), intent(in) :: unique_id
    character(str_def) :: axis_ref
    logical(l_def) :: has_axis
    call xios_is_defined_field_attr(unique_id, axis_ref=has_axis)
    if (has_axis) then
       call xios_get_field_attr(unique_id, axis_ref=axis_ref)
    else
       axis_ref = ""
    end if
  end function get_field_axis_ref

  !> @brief Return the dimension of an XIOS axis object.
  !> @param[in]    unique_id    XIOS id of the axis
  !> @return                    The dimension
  function get_axis_dimension(unique_id) result(dim)
    implicit none
    character(*), intent(in) :: unique_id
    integer(i_def) :: dim
    call xios_get_axis_attr(unique_id, n_glo=dim)
  end function get_axis_dimension

end module lfric_xios_diag_mod
