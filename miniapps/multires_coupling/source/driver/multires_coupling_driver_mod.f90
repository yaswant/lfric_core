!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the multires_coupling miniapp.
!>
!> This a base structure from which further physics-dynamics coupling development can
!> be done. A test algorithm can be called in order to test the scalar mapping and
!> setting up fields on different meshes. If the miniapp is not in test mode, the
!> regular gungho timestep is called on the dynamics mesh.
!>
module multires_coupling_driver_mod

  use multires_coupling_config_mod,             only : dynamics_mesh_name,          &
                                                       physics_mesh_name,           &
                                                       multires_coupling_mode,      &
                                                       multires_coupling_mode_test
  use variable_fields_mod,                      only : update_variable_fields
  use gungho_step_mod,                          only : gungho_step
  use clock_mod,                                only : clock_type
  use constants_mod,                            only : i_def, i_native, &
                                                       str_def, imdi
  use io_config_mod,                            only : write_diag,           &
                                                       use_xios_io,          &
                                                       diagnostic_frequency, &
                                                       nodal_output_on_w3
  use io_context_mod,                           only : io_context_type
  use mesh_collection_mod,                      only : mesh_collection
  use log_mod,                                  only : log_event,        &
                                                       LOG_LEVEL_ALWAYS, &
                                                       LOG_LEVEL_INFO
  use multires_coupling_mod,                    only : program_name
  use coupling_test_alg_mod,                    only : coupling_test_alg
  use xios,                                     only : xios_context_finalize
  use multires_coupling_model_mod,              only : initialise_model,          &
                                                       finalise_model,            &
                                                       initialise_infrastructure, &
                                                       finalise_infrastructure
  use gungho_model_data_mod,                    only : model_data_type,       &
                                                       create_model_data,     &
                                                       finalise_model_data,   &
                                                       initialise_model_data
  use multires_coupling_diagnostics_driver_mod, only : multires_coupling_diagnostics_driver
  use field_collection_mod,                     only : field_collection_type
  use field_mod,                                only : field_type

  implicit none

  private
  public initialise, run, finalise

  class(io_context_type), allocatable :: io_context

  ! Model run working data set
  type (model_data_type) :: dynamics_mesh_model_data
  type (model_data_type) :: physics_mesh_model_data

  ! Primary mesh ids
  integer(kind=i_def) :: prime_mesh_id = imdi
  integer(kind=i_def) :: prime_2D_mesh_id = imdi
  integer(kind=i_def) :: prime_shifted_mesh_id = imdi
  integer(kind=i_def) :: prime_double_level_mesh_id = imdi

  ! Dynamics and Physics mesh ids
  integer(kind=i_def) :: dynamics_mesh_id
  integer(kind=i_def) :: dynamics_2D_mesh_id
  integer(kind=i_def) :: physics_mesh_id
  integer(kind=i_def) :: physics_2D_mesh_id

  ! 2D Mesh names
  character(str_def) :: dynamics_2D_mesh_name
  character(str_def) :: physics_2D_mesh_name

  character(*), public, parameter   :: xios_ctx = 'multires_coupling'

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !>
  subroutine initialise( filename, model_communicator )

    implicit none

    character(:),      intent(in), allocatable :: filename
    integer(i_native), intent(in)              :: model_communicator

    class(clock_type), pointer :: clock

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------

    call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    !-------------------------------------------------------------------------
    ! Initialise Infrastructure such as meshes, FEM and runtime constants
    !-------------------------------------------------------------------------
    call log_event( 'Initialising Infrastructure...', LOG_LEVEL_ALWAYS )

    call initialise_infrastructure( model_communicator,        &
                                    filename,                  &
                                    program_name,              &
                                    io_context,                &
                                    prime_mesh_id,             &
                                    prime_2D_mesh_id,          &
                                    prime_shifted_mesh_id,     &
                                    prime_double_level_mesh_id )

    dynamics_2D_mesh_name = trim(dynamics_mesh_name)//'_2d'
    dynamics_mesh_id = mesh_collection%get_mesh_id(dynamics_mesh_name)
    dynamics_2D_mesh_id = mesh_collection%get_mesh_id(dynamics_2D_mesh_name)

    physics_2D_mesh_name = trim(physics_mesh_name)//'_2d'
    physics_mesh_id = mesh_collection%get_mesh_id(physics_mesh_name)
    physics_2D_mesh_id = mesh_collection%get_mesh_id(physics_2D_mesh_name)

    !-------------------------------------------------------------------------
    ! Create and initialise model data
    !-------------------------------------------------------------------------
    call log_event( 'Creating and Initialising Model Data...', LOG_LEVEL_ALWAYS )

    clock => io_context%get_clock()

    ! Instantiate the fields stored in model_data
    call create_model_data( dynamics_mesh_model_data, &
                            dynamics_mesh_id,         &
                            dynamics_2D_mesh_id,      &
                            clock )
    call create_model_data( physics_mesh_model_data, &
                            physics_mesh_id,         &
                            physics_2D_mesh_id,      &
                            clock )


    ! Initialise the fields stored in the model_data
    call initialise_model_data( dynamics_mesh_model_data, clock )
    call initialise_model_data( physics_mesh_model_data, clock )

   ! Initial output
   ! We only want these once at the beginning of a run
   if ( clock%is_initialisation() .and. write_diag .and. &
        multires_coupling_mode /= multires_coupling_mode_test ) then
     call multires_coupling_diagnostics_driver( dynamics_mesh_id,         &
                                                dynamics_2D_mesh_id,      &
                                                dynamics_mesh_model_data, &
                                                clock, nodal_output_on_w3 )
   end if

   if ( multires_coupling_mode /= multires_coupling_mode_test ) then
     ! Model configuration initialisation
     call initialise_model( clock,                    &
                            dynamics_mesh_id,         &
                            dynamics_mesh_model_data, &
                            physics_mesh_id,          &
                            physics_mesh_model_data )
   end if

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs time steps.
  !>
  subroutine run()

    implicit none

    class(clock_type), pointer :: clock

    clock => io_context%get_clock()

    if ( multires_coupling_mode == multires_coupling_mode_test ) then
      ! Call coupling test algorithm
      call coupling_test_alg()
    else
      ! Do timestep
      do while (clock%tick())
        ! Update time-varying fields
        call update_variable_fields( dynamics_mesh_model_data%ancil_times_list, &
                                      clock, dynamics_mesh_model_data%ancil_fields )
        ! Perform gungho timestep
        call gungho_step( dynamics_mesh_id,         &
                          dynamics_2D_mesh_id,      &
                          dynamics_mesh_model_data, &
                          clock )
        ! Write out output file
        call log_event(program_name//": Writing depository output", LOG_LEVEL_INFO)

        if ( (mod(clock%get_step(), diagnostic_frequency) == 0) &
              .and. (write_diag) ) then
              call multires_coupling_diagnostics_driver( dynamics_mesh_id,         &
                                                         dynamics_2D_mesh_id,      &
                                                         dynamics_mesh_model_data, &
                                                         clock, nodal_output_on_w3 )
        end if
      end do
    end if

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Finalise after run
  !>
  subroutine finalise()

    implicit none

    class(clock_type), pointer :: clock

    clock => io_context%get_clock()

    !-----------------------------------------------------------------------------
    ! Model finalise
    !-----------------------------------------------------------------------------
    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    !-------------------------------------------------------------------------
    ! Driver layer finalise
    !-------------------------------------------------------------------------

   ! Finalise XIOS context if we used it for diagnostic output or checkpointing
    if ( use_xios_io ) then
      call xios_context_finalize()
    end if

    if ( multires_coupling_mode /= multires_coupling_mode_test ) then
      call finalise_model( dynamics_mesh_id,         &
                           dynamics_mesh_model_data, &
                           physics_mesh_id,          &
                           physics_mesh_model_data,  &
                           program_name )
    end if

    ! Destroy the fields stored in model_data
    call finalise_model_data( dynamics_mesh_model_data )
    call finalise_model_data( physics_mesh_model_data )

    call finalise_infrastructure( program_name )

    call log_event( program_name//': Miniapp completed', LOG_LEVEL_INFO )

  end subroutine finalise

end module multires_coupling_driver_mod
