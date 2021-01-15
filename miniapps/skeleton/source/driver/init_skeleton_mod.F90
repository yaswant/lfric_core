!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Initialisation functionality for the skeleton miniapp

!> @details Handles init of prognostic fields and through the call to
!>          runtime_contants the coordinate fields and fem operators

module init_skeleton_mod

  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type
  use field_parent_mod,               only : write_interface
  use finite_element_config_mod,      only : element_order
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W3
  use log_mod,                        only : log_event,      &
                                             LOG_LEVEL_INFO, &
                                             LOG_LEVEL_ERROR
  use runtime_constants_mod,          only : create_runtime_constants
  use io_config_mod,                  only : write_diag, &
                                             use_xios_io
  use lfric_xios_write_mod,           only : write_field_face
  implicit none


  contains

  subroutine init_skeleton(mesh_id, twod_mesh_id, chi_xyz, chi_sph, panel_id, field_1)

    implicit none

    integer(i_def), intent(in)               :: mesh_id
    integer(i_def), intent(in)               :: twod_mesh_id
    ! Prognostic fields
    type( field_type ), intent(inout)        :: field_1
    ! Coordinate field
    type( field_type ), intent(inout)        :: chi_xyz(:)
    type( field_type ), intent(inout)        :: chi_sph(:)
    type( field_type ), intent(inout)        :: panel_id

    procedure(write_interface), pointer :: tmp_ptr

    call log_event( 'skeleton: Initialising miniapp ...', LOG_LEVEL_INFO )


    ! Create prognostic fields
    ! Creates a field in the W3 function space (fully discontinuous field)
    call field_1%initialise( vector_space = &
                      function_space_collection%get_fs(mesh_id, element_order, W3) )

    ! Set up field with an IO behaviour (XIOS only at present)

    if (write_diag .and. use_xios_io) then

       tmp_ptr => write_field_face

       call field_1%set_write_behaviour(tmp_ptr)

    end if

    ! Create runtime_constants object. This in turn creates various things
    ! needed by the fem algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_runtime_constants(mesh_id, twod_mesh_id, chi_xyz, chi_sph, panel_id)

    call log_event( 'skeleton: Miniapp initialised', LOG_LEVEL_INFO )

  end subroutine init_skeleton

end module init_skeleton_mod
