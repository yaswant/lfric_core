!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page Miniapp da_dev program

!> @brief Main program used to illustrate how to write LFRic miniapps.

!> @details Calls init, run and finalise routines from a driver module

program da_dev

  use da_dev_driver_mod, only : initialise, run, finalise

  implicit none

  call initialise()

  call run()

  call finalise()

end program da_dev
