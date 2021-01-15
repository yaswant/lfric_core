!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Create the prognostic fields and place them in the depository

!> @details Create the prognostic fields and place them both in the
!>          depository field collection and put pointers to them in the
!>          prognostic field collection

module create_gravity_wave_prognostics_mod

  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type
  use field_parent_mod,               only : write_interface, &
                                             checkpoint_write_interface, &
                                             checkpoint_read_interface
  use field_collection_mod,           only : field_collection_type
  use finite_element_config_mod,      only : element_order
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W0, W2, W3, Wtheta
  use io_config_mod,                  only : write_diag, &
                                             use_xios_io, &
                                             checkpoint_write, &
                                             checkpoint_read
  use io_mod,                         only : checkpoint_read_netcdf, &
                                             checkpoint_write_netcdf
  use lfric_xios_read_mod,            only : checkpoint_read_xios
  use lfric_xios_write_mod,           only : write_field_node, &
                                             write_field_face, &
                                             checkpoint_write_xios

  use gravity_wave_constants_config_mod,&
                                      only : b_space, &
                                             b_space_w0, &
                                             b_space_w3, &
                                             b_space_wtheta
  use log_mod,                        only : log_event, &
                                             LOG_LEVEL_INFO, &
                                             LOG_LEVEL_ERROR
  use pure_abstract_field_mod,        only : pure_abstract_field_type
  implicit none

  private
  public create_gravity_wave_prognostics

  contains

  !> @brief Create the prognostic fields and place them in the depository
  !> @param [in] mesh_id The identifier of the primary mesh
  !> @param [inout] depository A collection of all fields used in the miniapp
  !> @param [inout] prognostics A collection of all the prognostic fields.
  !>                            For the gravity wave miniapp, the prognostics
  !>                            are all the fields in the depository - so the
  !>                            two collections contain the same fields
  subroutine create_gravity_wave_prognostics(mesh_id, depository, prognostics)

    implicit none
    integer(i_def), intent(in)                   :: mesh_id
    type( field_collection_type ), intent(inout) :: depository
    type( field_collection_type ), intent(inout) :: prognostics

    type( field_type) :: wind
    type( field_type) :: buoyancy
    type( field_type) :: pressure

    integer(i_def) :: buoyancy_space

    class(pure_abstract_field_type), pointer  :: tmp_ptr => null()

    procedure(write_interface),       pointer :: tmp_write_ptr
    procedure(checkpoint_write_interface), pointer :: tmp_checkpoint_write_ptr
    procedure(checkpoint_read_interface),  pointer :: tmp_checkpoint_read_ptr

    ! Create prognostic fields
    select case(b_space)
      case(b_space_w0)
        buoyancy_space = W0
        call log_event( 'gravity_wave: Using W0 for buoyancy', LOG_LEVEL_INFO )
      case(b_space_w3)
        buoyancy_space = W3
        call log_event( 'gravity_wave: Using W3 for buoyancy', LOG_LEVEL_INFO )
      case(b_space_wtheta)
        buoyancy_space = Wtheta
        call log_event( 'gravity_wave: Using Wtheta for buoyancy', LOG_LEVEL_INFO )
      case default
        call log_event( 'gravity_wave: Invalid buoyancy space', LOG_LEVEL_ERROR )
    end select

    call wind%initialise( vector_space = &
                       function_space_collection%get_fs(mesh_id, element_order, W2), &
                       name="wind" )
    call buoyancy%initialise( vector_space = &
               function_space_collection%get_fs(mesh_id, element_order, buoyancy_space), &
               name="buoyancy" )
    call pressure%initialise( vector_space = &
                       function_space_collection%get_fs(mesh_id, element_order, W3), &
                       name="pressure" )


    ! Set I/O behaviours for diagnostic output
    if (write_diag .and. use_xios_io) then
       ! Fields that are output on the XIOS face domain
       tmp_write_ptr => write_field_face
       call wind%set_write_behaviour(tmp_write_ptr)
       call pressure%set_write_behaviour(tmp_write_ptr)
       if (buoyancy_space == W0) then
         tmp_write_ptr => write_field_node
         call buoyancy%set_write_behaviour(tmp_write_ptr)
       else
         tmp_write_ptr => write_field_face
         call buoyancy%set_write_behaviour(tmp_write_ptr)
       end if
    end if

    ! Set I/O behaviours for checkpoint / restart
    if ( checkpoint_write .or. checkpoint_read) then
      if (use_xios_io) then
        ! Use XIOS for checkpoint / restart
        tmp_checkpoint_write_ptr => checkpoint_write_xios
        tmp_checkpoint_read_ptr => checkpoint_read_xios
        call log_event( 'GungHo: Using XIOS for checkpointing...', LOG_LEVEL_INFO )
      else
        ! Use old checkpoint and restart methods
        tmp_checkpoint_write_ptr => checkpoint_write_netcdf
        tmp_checkpoint_read_ptr => checkpoint_read_netcdf
        call log_event( 'GungHo: Using NetCDF for checkpointing...', LOG_LEVEL_INFO )
      end if

      call wind%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)
      call pressure%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)
      call buoyancy%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)

      call wind%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)
      call pressure%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)
      call buoyancy%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)

    end if

    ! Put the prognostic fields into the depository
    call depository%add_field(wind)
    call depository%add_field(pressure)
    call depository%add_field(buoyancy)

    ! Put pointers to the prognostic fields into the prognostic collection
    tmp_ptr => depository%get_field('wind')
    call prognostics%add_reference_to_field(tmp_ptr)
    tmp_ptr => depository%get_field('pressure')
    call prognostics%add_reference_to_field(tmp_ptr)
    tmp_ptr => depository%get_field('buoyancy')
    call prognostics%add_reference_to_field(tmp_ptr)

  end subroutine create_gravity_wave_prognostics

end module create_gravity_wave_prognostics_mod
