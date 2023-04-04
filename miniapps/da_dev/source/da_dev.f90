!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page Miniapp da_dev program

!> @brief Main program for running da_dev independently.

!> @details Calls init, run and finalise routines from a driver module

program da_dev

  use da_dev_driver_mod,         only : initialise_lfric, run, finalise_lfric, &
                                        initialise_model, finalise_model,      &
                                        initialise_lfric_comm
  use driver_model_data_mod,     only : model_data_type
  use constants_mod,             only : i_native

  implicit none

  type(model_data_type) :: model_data
  integer(i_native)     :: model_communicator

  character(*), parameter :: program_name = "da_dev"

  call initialise_lfric_comm(program_name, model_communicator)

  call initialise_lfric(program_name, model_communicator)

  call initialise_model(program_name, model_communicator, model_data)

  call run(program_name, model_data)

  call finalise_model(program_name, model_data%depository)

  call finalise_lfric(program_name)

end program da_dev
