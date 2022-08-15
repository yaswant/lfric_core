!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the da_dev miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module da_dev_driver_mod

  use checksum_alg_mod,           only : checksum_alg
  use cli_mod,                    only : get_initial_filename
  use clock_mod,                  only : clock_type
  use configuration_mod,          only : final_configuration
  use constants_mod,              only : i_def, i_native, &
                                         PRECISION_REAL, r_def
  use convert_to_upper_mod,       only : convert_to_upper
  use driver_comm_mod,            only : init_comm, final_comm
  use driver_model_data_mod,      only : model_data_type
  use driver_log_mod,             only : init_logger, final_logger
  use driver_mesh_mod,            only : init_mesh, final_mesh
  use driver_fem_mod,             only : init_fem, final_fem
  use driver_io_mod,              only : init_io, final_io, &
                                         get_clock
  use field_mod,                  only : field_type
  use field_collection_mod,       only : field_collection_type
  use init_da_dev_mod,            only : init_da_dev
  use io_config_mod,              only : write_diag
  use log_mod,                    only : log_event,          &
                                         log_scratch_space,  &
                                         LOG_LEVEL_ALWAYS,   &
                                         LOG_LEVEL_INFO
  use mesh_mod,                   only : mesh_type
  use mpi_mod,                    only : get_comm_size, &
                                         get_comm_rank
  use da_dev_mod,                 only : load_configuration
  use da_dev_alg_mod,             only : da_dev_alg

  implicit none

  private
  public initialise, run, finalise

  character(*), parameter :: program_name = "da_dev"

  ! Prognostic fields
  type(model_data_type)                  :: model_data

  ! Coordinate field
  type(field_type), target, dimension(3) :: chi
  type(field_type), target               :: panel_id
  type(mesh_type),  pointer              :: mesh      => null()
  type(mesh_type),  pointer              :: twod_mesh => null()

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !>
  subroutine initialise()

    implicit none

    character(:), allocatable :: filename
    integer(i_native) :: model_communicator

    class(clock_type),      pointer :: clock => null()
    real(r_def)                     :: dt_model

    call init_comm( program_name, model_communicator )

    call get_initial_filename( filename )
    call load_configuration( filename, program_name )

    call init_logger(get_comm_size(), get_comm_rank(), program_name)

    write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------
    call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Create the mesh
    call init_mesh( get_comm_rank(), get_comm_size(), mesh, &
                    twod_mesh = twod_mesh )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh, chi, panel_id )

    ! Initialise I/O context
    call init_io( program_name, model_communicator, chi, panel_id )

    ! Get dt from model clock
    clock => get_clock()
    dt_model = real(clock%get_seconds_per_step(), r_def)

    ! Create and initialise prognostic fields
    call init_da_dev(mesh, chi, panel_id, dt_model, model_data)

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs time steps.
  !>
  subroutine run()

    implicit none

    class(clock_type), pointer :: clock => null()
    logical                    :: running
    type(field_type), pointer  :: working_field => null()

    clock => get_clock()
    running = clock%tick()
    working_field => model_data%depository%get_field( 'da_dev_field' )

    ! Call an algorithm
    call da_dev_alg( working_field )

    ! Write out output file
    call log_event(program_name//": Writing diagnostic output", LOG_LEVEL_INFO)

    if ( write_diag ) then
      ! Calculation and output of diagnostics
      call working_field%write_field('da_dev_field')
    end if

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tidies up after a run.
  !>
  subroutine finalise()

    implicit none

    type(field_type), pointer :: working_field => null()

    working_field => model_data%depository%get_field( 'da_dev_field' )

    !---------------------------------------------------------------------------
    ! Model finalise
    !---------------------------------------------------------------------------
    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Write checksums to file
    call checksum_alg(program_name, working_field, 'da_dev_field')

    call log_event( program_name//': Miniapp completed', LOG_LEVEL_INFO )

    !-------------------------------------------------------------------------
    ! Driver layer finalise
    !-------------------------------------------------------------------------

    ! Finalise IO
    call final_io()

    call final_fem()

    call final_mesh()

    call final_configuration()

    call final_logger( program_name )

    call final_comm()

  end subroutine finalise

end module da_dev_driver_mod
