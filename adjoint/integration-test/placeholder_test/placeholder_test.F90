!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief A placeholder integration test so that the app can function.
!> @details Adapted from https://code.metoffice.gov.uk/trac/lfric/wiki/LFRicTechnical/Testing/IntegrationTesting
!>          and linear integration tests.
program placeholder_test

  use configuration_mod,      only: read_configuration, final_configuration
  use constants_mod,          only: r_def
  use driver_collections_mod, only: init_collections, final_collections
  use driver_time_mod,        only: init_time, get_calendar
  use gungho_modeldb_mod,     only: modeldb_type
  use halo_comms_mod,         only: initialise_halo_comms, &
                                    finalise_halo_comms
  use log_mod,                only: log_event,       &
                                    LOG_LEVEL_ERROR, &
                                    LOG_LEVEL_INFO,  &
                                    log_scratch_space
  use mpi_mod,                only: create_comm, destroy_comm, global_mpi

  implicit none

  ! Model run working data set
  type(modeldb_type) :: modeldb

  character(*), parameter :: application_name = "placeholder_test"

  character(:), allocatable :: filename

  ! Variables used for parsing command line arguments
  integer :: length, status, nargs
  character(len=0) :: dummy
  character(len=:), allocatable :: program_name

  integer :: communicator

  ! Result of computation
  real(kind=r_def) :: res

  modeldb%mpi => global_mpi

  call create_comm( communicator )
  call modeldb%mpi%initialise( communicator )
  call initialise_halo_comms( communicator )

  call log_event( 'Placeholder test running ...', LOG_LEVEL_INFO )

  ! Create the depository, prognostics and diagnostics field collections
  call modeldb%fields%add_empty_field_collection("depository", table_len = 100)
  call modeldb%model_data%prognostic_fields%initialise(name="prognostics", table_len=100)
  call modeldb%model_data%diagnostic_fields%initialise(name="diagnostics", table_len=100)

  ! Parse command line parameters
  call get_command_argument( 0, dummy, length, status )
  allocate(character(length)::program_name)
  call get_command_argument( 0, program_name, length, status )
  nargs = command_argument_count()

  call get_command_argument( 1, dummy, length, status )
  allocate( character(length) :: filename )
  call get_command_argument( 1, filename, length, status )

  call modeldb%configuration%initialise( program_name, table_len=10 )
  call read_configuration( filename, modeldb%configuration )
  deallocate( filename )

  call init_collections()
  call init_time( modeldb%clock )

  ! Some science is called, resulting in the following result
  res = 0.0_r_def
  write( log_scratch_space, '("Result = ", F5.2)' ) res
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  call final_collections()
  call final_configuration()
  call finalise_halo_comms()
  call destroy_comm()

end program placeholder_test
