!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief   The top level program for the da dev integration tests.
!>@details Sets up and runs the integration tests specified in
!!         algorithms_test.py. Currently only
!!         da_dev_increment_alg_mod.x90 is tested, this algorithm
!!         takes a real field and adds one to its value.
program algorithms_test

  use constants_mod,                 only : i_def, r_def
  use halo_comms_mod,                only : initialise_halo_comms, &
                                            finalise_halo_comms
  use log_mod,                       only : log_event,          &
                                            initialise_logging, &
                                            finalise_logging,   &
                                            LOG_LEVEL_ERROR,    &
                                            LOG_LEVEL_INFO
  use mpi_mod,                       only : initialise_comm, store_comm, &
                                            finalise_comm,               &
                                            get_comm_size, get_comm_rank

  ! For the da_dev_increment_alg_mod test
  use da_dev_increment_alg_mod,      only : da_dev_increment_alg
  use da_dev_mod,                    only : load_configuration
  use driver_mesh_mod,               only : init_mesh, final_mesh
  use driver_fem_mod,                only : init_fem, final_fem
  use field_mod,                     only : field_type, field_proxy_type
  use fs_continuity_mod,             only : Wtheta
  use function_space_collection_mod, only : function_space_collection
  use function_space_mod,            only : function_space_type
  use mesh_mod,                      only : mesh_type

  implicit none

  ! MPI communicator
  integer(i_def) :: comm

  ! Number of processes and local rank
  integer(i_def) :: total_ranks, local_rank

  character(:), allocatable :: filename

  ! Variables used for parsing command line arguments
  integer :: length, status, nargs
  character(len=0) :: dummy
  character(len=:), allocatable :: program_name, test_flag

  ! Flags which determine the tests that will be carried out
  logical :: do_test_da_dev_increment_alg_mod = .false.

  ! Usage message to print
  character(len=256) :: usage_message

  ! For the da_dev_increment_alg_mod test
  ! Coordinate field
  type(field_type), target, dimension(3) :: chi
  type(field_type), target               :: panel_id
  type(mesh_type),  pointer              :: mesh      => null()
  type(function_space_type), pointer     :: wtheta_fs => null()
  integer(i_def),   parameter            :: element_order = 0
  ! Test Field
  type(field_type), target               :: test_field
  type(field_proxy_type)                 :: field_proxy

  character(len=128)                     :: output
  real(r_def), parameter                 :: tol = 1.0e-3_r_def

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Communicators and Logging Setup
  ! Initialise MPI communicatios and get a valid communicator
  call initialise_comm(comm)

  ! Save the communicator for later use
  call store_comm(comm)

  ! Initialise halo functionality
  call initialise_halo_comms(comm)

  total_ranks = get_comm_size()
  local_rank  = get_comm_rank()

  call initialise_logging( comm, 'da_dev_test' )

  call log_event( 'da dev testing running ...', LOG_LEVEL_INFO )

  ! Parse command line parameters
  call get_command_argument( 0, dummy, length, status )
  allocate(character(length)::program_name)
  call get_command_argument( 0, program_name, length, status )
  nargs = command_argument_count()

  ! Print out usage message if wrong number of arguments is specified
  if (nargs /= 2) then
    write(usage_message,*) "Usage: ",trim(program_name), &
      " <namelist filename> "       // &
      " test_XXX with XXX in { "    // &
      " da_dev_increment_alg_mod, " // &
      " } "
    call log_event( trim(usage_message), LOG_LEVEL_ERROR )
  end if

  call get_command_argument( 1, dummy, length, status )
  allocate( character(length) :: filename )
  call get_command_argument( 1, filename, length, status )

  call get_command_argument( 2, dummy, length, status )
  allocate(character(length)::test_flag)
  call get_command_argument( 2, test_flag, length, status )

  ! Choose test case depending on flag provided in the first command
  ! line argument
  select case (trim(test_flag))
  case ("test_da_dev_increment_alg_mod")
    do_test_da_dev_increment_alg_mod = .true.
  case default
    call log_event( "Unknown test", LOG_LEVEL_ERROR )
  end select

  if (do_test_da_dev_increment_alg_mod) then
    ! Load config and create mesh and function space
    call load_configuration( filename, program_name )
    call init_mesh( local_rank, total_ranks, mesh )
    call init_fem( mesh, chi, panel_id )
    wtheta_fs => function_space_collection%get_fs( mesh,          &
                                                   element_order, &
                                                   Wtheta )
    ! Initialise test field and get proxy
    call test_field%initialise( wtheta_fs, name='test_field' )
    field_proxy = test_field%get_proxy()
    ! Set first dof to nominal value
    field_proxy%data(:) = 1.0_r_def

    ! Call algorithm
    call da_dev_increment_alg(test_field)

    write ( output, '("Returned Field Values =", ES15.8)' ) field_proxy%data(1)

    ! Test if value is correct
    if ( all( abs(field_proxy%data(:) - 2.0_r_def)  <= tol ) ) then
      call log_event( 'test PASS', LOG_LEVEL_INFO )
      call log_event( output, LOG_LEVEL_INFO )
    else
      call log_event( 'test FAIL', LOG_LEVEL_INFO )
      call log_event( 'Expected Field Values = 2.00000000E+00', LOG_LEVEL_INFO )
      call log_event( output, LOG_LEVEL_INFO )
    end if

    ! Finalise
    call final_fem()
    call final_mesh()
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Finalise and close down

  if (allocated(program_name)) deallocate(program_name)
  if (allocated(filename))     deallocate(filename)
  if (allocated(test_flag))    deallocate(test_flag)
  nullify(mesh)
  nullify(wtheta_fs)

  call log_event( 'da dev functional testing completed ...', LOG_LEVEL_INFO )

  ! Finalise halo functionality
  call finalise_halo_comms()

  ! Finalise MPI communications
  call finalise_comm()

  ! Finalise the logging system
  call finalise_logging()

end program algorithms_test
