!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Handles calling of calender update for gungho

module gungho_update_calendar_mod

  use time_config_mod, only : timestep_start
  use log_mod,         only : log_event, &
                              LOG_LEVEL_INFO
  use io_config_mod,   only : use_xios_io
  use xios,            only : xios_update_calendar
  use constants_mod,   only : i_def

  implicit none

  private
  public gungho_update_calendar

  contains

  !> @brief Call calendar update if required
  !> @param[in] timestep number of current timestep
  subroutine gungho_update_calendar( timestep )

    implicit none

    integer(i_def), intent(in) :: timestep

    ! Update XIOS calendar if we are using it for diagnostic output or
    ! checkpoint
    if ( use_xios_io ) then
      call log_event( "Gungho: Updating XIOS timestep", LOG_LEVEL_INFO )
      call xios_update_calendar(timestep - timestep_start + 1)
    end if

  end subroutine gungho_update_calendar

end module gungho_update_calendar_mod
