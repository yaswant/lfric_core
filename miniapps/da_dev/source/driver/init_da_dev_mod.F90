!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Initialisation functionality for the da_dev miniapp

!> @details Handles init of prognostic fields and through the call to
!>          runtime_contants the coordinate fields and fem operators

module init_da_dev_mod

  use constants_mod,                  only : i_def, r_def
  use driver_model_data_mod,          only : model_data_type
  use field_mod,                      only : field_type
  use field_collection_mod,           only : field_collection_type
  use field_parent_mod,               only : write_interface
  use finite_element_config_mod,      only : element_order
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W3
  use log_mod,                        only : log_event,      &
                                             LOG_LEVEL_INFO, &
                                             LOG_LEVEL_ERROR
  use mesh_mod,                       only : mesh_type
  use pure_abstract_field_mod,        only : pure_abstract_field_type
  use io_config_mod,                  only : write_diag, &
                                             use_xios_io
  use lfric_xios_write_mod,           only : write_field_face
  use da_dev_constants_mod,           only : create_da_dev_constants

  implicit none

  contains

  subroutine init_da_dev( mesh, chi, panel_id, dt, model_data)

    implicit none

    type(mesh_type), intent(in), pointer        :: mesh
    real(r_def),    intent(in)                  :: dt

    ! Prognostic fields
    type( model_data_type ), intent(inout)      :: model_data
    type( field_type ), allocatable, target     :: field

    ! Coordinate field
    type( field_type ), intent(inout)           :: chi(:)
    type( field_type ), intent(inout)           :: panel_id

    procedure(write_interface), pointer         :: tmp_ptr

    call log_event( 'da_dev: Initialising miniapp ...', LOG_LEVEL_INFO )

    call model_data%depository%initialise(name = 'depository', table_len=100)

    allocate(field)

    ! Create prognostic fields
    ! Creates a field in the W3 function space (fully discontinuous field)
    call field%initialise( &
      vector_space = function_space_collection%get_fs(mesh, element_order, W3), &
      name = 'da_dev_field' )

    ! Set up field with an IO behaviour (XIOS only at present)
    if (write_diag .and. use_xios_io) then
       tmp_ptr => write_field_face
       call field%set_write_behaviour(tmp_ptr)

    end if

    call model_data%depository%add_field(field)

    ! Create da_dev runtime constants. This creates various things
    ! needed by the fem algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_da_dev_constants(mesh, chi, panel_id)

    call log_event( 'da_dev: Miniapp initialised', LOG_LEVEL_INFO )

  end subroutine init_da_dev

end module init_da_dev_mod
