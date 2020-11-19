!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief init functionality for demonstrating catalyst

!> @details Handles init of prognostic and coordinate fields

module init_catalyst_demo_mod

  use assign_coordinate_field_mod,    only : assign_coordinate_field
  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type, &
                                             write_diag_interface
  use field_parent_mod,               only : checkpoint_write_interface, &
                                             checkpoint_read_interface
  use finite_element_config_mod,      only : element_order

  use fs_continuity_mod,              only : W0, W2, W3, Wtheta
  use gw_init_fields_alg_mod,         only : gw_init_fields_alg
  use log_mod,                        only : log_event, log_scratch_space, &
                                             LOG_LEVEL_INFO, &
                                             LOG_LEVEL_ERROR
  use gravity_wave_constants_config_mod,only : b_space,                           &
                                               gravity_wave_constants_b_space_w0, &
                                               gravity_wave_constants_b_space_w3, &
                                               gravity_wave_constants_b_space_wtheta
  use runtime_constants_mod,          only : create_runtime_constants
  use io_config_mod,                  only : write_diag,       &
                                             use_xios_io,      &
                                             checkpoint_write, &
                                             checkpoint_read
  use write_methods_mod,              only : write_field_face, &
                                             write_field_node, &
                                             checkpoint_write_xios, &
                                             checkpoint_write_netcdf
  use read_methods_mod,               only : checkpoint_read_netcdf, &
                                             checkpoint_read_xios

  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection
  use function_space_chain_mod,       only : function_space_chain_type

  use formulation_config_mod,         only : l_multigrid
  use multigrid_config_mod,           only : multigrid_chain_nitems

  implicit none


  contains

  subroutine init_catalyst_demo( mesh_id, twod_mesh_id, multigrid_mesh_ids, &
                                 chi, multigrid_function_space_chain,       &
                                 wind, pressure, buoyancy )

    implicit none

    integer(i_def),  intent(in)  :: mesh_id
    integer(i_def),  intent(in)  :: twod_mesh_id
    integer(i_def),  intent(in)  :: multigrid_mesh_ids(:)

    type(function_space_chain_type), intent(out) :: multigrid_function_space_chain

    ! Prognostic fields
    type( field_type ), intent(inout)        :: wind, pressure, buoyancy
    ! Coordinate field
    type( field_type ), intent(inout)        :: chi(:)
    integer(i_def)                           :: buoyancy_space

    integer(i_def) :: i

    procedure(write_diag_interface), pointer :: tmp_write_diag_ptr
    procedure(checkpoint_write_interface), pointer :: tmp_checkpoint_write_ptr
    procedure(checkpoint_read_interface), pointer  :: tmp_checkpoint_read_ptr

    type(function_space_type),  pointer :: function_space => null()

    call log_event( 'catalyst_demo: Initialising miniapp ...', LOG_LEVEL_INFO )

    !===============================================================================
    ! Now create the function space chain for multigrid
    !===============================================================================
    if (l_multigrid) then

      do i=1, size(multigrid_mesh_ids)

        ! Make sure this function_space is in the collection
        function_space => function_space_collection%get_fs( multigrid_mesh_ids(i), &
                                                            0,                     &
                                                            W3 )

        write( log_scratch_space,"(A,I0,A)")                       &
             'Adding function_space id ', function_space%get_id(), &
             ' to multigrid function_space chain'
        call log_event( log_scratch_space, LOG_LEVEL_INFO )

        call multigrid_function_space_chain%add( function_space )

      end do
    end if


    ! Create prognostic fields
    select case(b_space)
      case(gravity_wave_constants_b_space_w0)
        buoyancy_space = W0
        call log_event( 'catalyst_demo: Using W0 for buoyancy', LOG_LEVEL_INFO )
      case(gravity_wave_constants_b_space_w3)
        buoyancy_space = W3
        call log_event( 'catalyst_demo: Using W3 for buoyancy', LOG_LEVEL_INFO )
      case(gravity_wave_constants_b_space_wtheta)
        buoyancy_space = Wtheta
        call log_event( 'catalyst_demo: Using Wtheta for buoyancy', LOG_LEVEL_INFO )
      case default
        call log_event( 'catalyst_demo: Invalid buoyancy space', LOG_LEVEL_ERROR )
    end select

    call wind%initialise( vector_space = &
                       function_space_collection%get_fs(mesh_id, element_order, W2))
    call buoyancy%initialise( vector_space = &
               function_space_collection%get_fs(mesh_id, element_order, buoyancy_space) )
    call pressure%initialise( vector_space = &
                       function_space_collection%get_fs(mesh_id, element_order, W3) )

    ! Set I/O behaviours for diagnostic output

    if (write_diag .and. use_xios_io) then

       ! Fields that are output on the XIOS face domain

       tmp_write_diag_ptr => write_field_face

       call wind%set_write_diag_behaviour(tmp_write_diag_ptr)
       call pressure%set_write_diag_behaviour(tmp_write_diag_ptr)

       if (buoyancy_space == W0) then
         tmp_write_diag_ptr => write_field_node
         call buoyancy%set_write_diag_behaviour(tmp_write_diag_ptr)
       else
         tmp_write_diag_ptr => write_field_face
         call buoyancy%set_write_diag_behaviour(tmp_write_diag_ptr)
       end if

    end if

    ! Set I/O behaviours for checkpoint / restart

    if ( checkpoint_write .or. checkpoint_read) then

      if (use_xios_io) then

        ! Use XIOS for checkpoint / restart

        tmp_checkpoint_write_ptr => checkpoint_write_xios
        tmp_checkpoint_read_ptr => checkpoint_read_xios

        call log_event( 'catalyst_demo: Using XIOS for checkpointing...', LOG_LEVEL_INFO )

      else

        ! Use old checkpoint and restart methods

        tmp_checkpoint_write_ptr => checkpoint_write_netcdf
        tmp_checkpoint_read_ptr => checkpoint_read_netcdf

        call log_event( 'catalyst_demo: Using NetCDF for checkpointing...', LOG_LEVEL_INFO )

      end if

      call wind%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)
      call pressure%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)
      call buoyancy%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)

      call wind%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)
      call pressure%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)
      call buoyancy%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)

    end if

    ! Create runtime_constants object. This in turn creates various things
    ! needed by the timestepping algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_runtime_constants(mesh_id, twod_mesh_id, chi)

    ! Initialise prognostic fields
    call gw_init_fields_alg(wind, pressure, buoyancy)

    nullify( tmp_write_diag_ptr, tmp_checkpoint_write_ptr, &
             tmp_checkpoint_read_ptr, function_space )

    call log_event( 'catalyst_demo: Miniapp initialised', LOG_LEVEL_INFO )

  end subroutine init_catalyst_demo

end module init_catalyst_demo_mod
