!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>@brief Drives the execution of the (tangent) linear model.
module linear_driver_mod

  use clock_mod,                  only : clock_type
  use constants_mod,              only : i_def, i_native, imdi
  use field_mod,                  only : field_type
  use gungho_mod,                 only : program_name
  use gungho_diagnostics_driver_mod, &
                                  only : gungho_diagnostics_driver
  use gungho_model_mod,           only : initialise_infrastructure, &
                                         initialise_model, &
                                         finalise_infrastructure, &
                                         finalise_model
  use gungho_model_data_mod,      only : model_data_type, &
                                         create_model_data, &
                                         initialise_model_data, &
                                         output_model_data, &
                                         finalise_model_data
  use gungho_step_mod,            only : gungho_step
  use io_config_mod,              only : write_diag, &
                                         diagnostic_frequency, &
                                         nodal_output_on_w3
  use io_context_mod,             only : io_context_type
  use log_mod,                    only : log_event,         &
                                         log_scratch_space, &
                                         LOG_LEVEL_ALWAYS
  use linear_model_data_mod,      only : linear_create_ls, &
                                         linear_init_ls,   &
                                         linear_init_pert
  use linear_model_mod,           only : initialise_linear_model, &
                                         finalise_linear_model
  use linear_diagnostics_driver_mod, &
                                  only : linear_diagnostics_driver
  use linear_step_mod,            only : linear_step
  use linear_data_algorithm_mod,  only : update_ls_file_alg
  use initialization_config_mod,  only : ls_option, &
                                         ls_option_file

  implicit none

  private
  public initialise, run, finalise

  ! Model run working data set
  type (model_data_type) :: model_data

  integer(i_def) :: mesh_id              = imdi
  integer(i_def) :: twod_mesh_id         = imdi
  integer(i_def) :: shifted_mesh_id      = imdi
  integer(i_def) :: double_level_mesh_id = imdi

  class(io_context_type), allocatable :: io_context

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Sets up required state in preparation for run.
  !>@param[in] filename Name of the file containing the desired configuration
  !>@param[in] model_communicator The MPI communicator for use within the model
  subroutine initialise( filename, model_communicator )

    implicit none

    character(*),      intent(in) :: filename
    integer(i_native), intent(in) :: model_communicator

    class(clock_type), pointer :: clock

    ! Initialise infrastructure and setup constants
    call initialise_infrastructure( model_communicator,   &
                                    filename,             &
                                    program_name,         &
                                    io_context,           &
                                    mesh_id,              &
                                    twod_mesh_id,         &
                                    shifted_mesh_id,      &
                                    double_level_mesh_id, &
                                    model_data )

    clock => io_context%get_clock()
    ! Instantiate the fields stored in model_data
    call create_model_data( model_data,   &
                            mesh_id,      &
                            twod_mesh_id, &
                            clock )

    ! Instantiate the linearisation state
    call linear_create_ls( model_data,   &
                           mesh_id,      &
                           twod_mesh_id )

    ! Initialise the fields stored in the model_data
    call initialise_model_data( model_data, clock )

    ! Model configuration initialisation
    call initialise_model( clock,   &
                           mesh_id, &
                           model_data )

    ! Initialise the linearisation state
    call linear_init_ls( mesh_id,      &
                         twod_mesh_id, &
                         model_data,   &
                         clock )

    ! Initialise the linear model perturbation state
    call linear_init_pert( mesh_id,      &
                           twod_mesh_id, &
                           model_data )

    ! Linear model configuration initialisation
    call initialise_linear_model( clock,   &
                                  mesh_id, &
                                  model_data )

    ! Initial output
    if ( clock%is_initialisation() .and. write_diag ) then
        ! Calculation and output of diagnostics
        call gungho_diagnostics_driver( mesh_id,      &
                                        twod_mesh_id, &
                                        model_data,   &
                                        clock,        &
                                        nodal_output_on_w3 )

        call linear_diagnostics_driver( mesh_id,    &
                                        model_data, &
                                        clock,      &
                                        nodal_output_on_w3 )
    end if



  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Timesteps the model, calling the desired timestepping algorithm based
  !upon the configuration
  subroutine run()

    implicit none

    class(clock_type), pointer :: clock

    write( log_scratch_space,'(A)' ) 'Running '//program_name//' ...'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    clock => io_context%get_clock()
    do while ( clock%tick() )

      if ( ls_option == ls_option_file ) then
        call update_ls_file_alg( model_data%ls_times_list, &
                                 clock,                    &
                                 model_data%ls_fields,     &
                                 model_data%ls_mr,         &
                                 model_data%ls_moist_dyn )
      end if

      call linear_step( mesh_id,      &
                        twod_mesh_id, &
                        model_data,   &
                        clock )

      if ( ( mod(clock%get_step(), diagnostic_frequency) == 0 ) &
           .and. ( write_diag ) ) then

        ! Calculation and output diagnostics
        call gungho_diagnostics_driver( mesh_id,      &
                                        twod_mesh_id, &
                                        model_data,   &
                                        clock,        &
                                        nodal_output_on_w3 )

        call linear_diagnostics_driver( mesh_id,    &
                                        model_data, &
                                        clock,      &
                                        nodal_output_on_w3 )
      end if

    end do ! end ts loop

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tidies up after a run.
  subroutine finalise()

    implicit none

    class(clock_type), pointer :: clock

    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    clock => io_context%get_clock()

    ! Output the fields stored in the model_data (checkpoint and dump)
    call output_model_data( model_data, clock )

    ! Model configuration finalisation
    call finalise_model( mesh_id,    &
                         model_data, &
                         program_name )

    call finalise_linear_model( )

    ! Destroy the fields stored in model_data
    call finalise_model_data( model_data )

    ! Finalise infrastructure and constants
    call finalise_infrastructure( program_name )

  end subroutine finalise

end module linear_driver_mod
