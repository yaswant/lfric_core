!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Functions to compute a variety of analytic profiles
!> @details Collection of functions to return the value of a field at a given
!!          point based upon a specified analytic formula
module analytic_streamfunction_profiles_mod

use constants_mod,           only: r_def, i_def, pi
use initial_wind_config_mod, only: profile_sbr_streamfunction,    &
                                   profile_dcmip301_streamfunction
use planet_config_mod,       only: scaled_radius
use log_mod,                 only: log_event,                &
                                   log_scratch_space,        &
                                   LOG_LEVEL_ERROR

implicit none

private
public :: analytic_streamfunction

contains

!> @brief Compute an analytic streamfunction field
!> @param[in] chi Position in physical coordinates
!> @param[in] choice Integer defining which specified formula to use
!> @param[in] num_options Number of scalar options to supply
!> @param[in] option Array of real values used to generate the initial profile
!> @return psi Result streamfunction field vector u = curl(psi)
function analytic_streamfunction(chi, choice, num_options, option) result(psi)

  implicit none

  real(kind=r_def),    intent(in) :: chi(3)
  integer(kind=i_def), intent(in) :: choice, num_options
  real(kind=r_def),    optional   :: option(num_options)
  real(kind=r_def)                :: psi(3)
  real(kind=r_def)                :: s
  real(kind=r_def)                :: lat_pole, lon_pole

  if ( .not. present(option) ) option(:) = 0.0_r_def

  select case ( choice )

    case ( profile_sbr_streamfunction, & 
           profile_dcmip301_streamfunction)      
      s = 0.5_r_def*(chi(3)/scaled_radius + 1.0_r_def)
      ! Turn off the height variation for the dcmip test
      if ( choice == profile_dcmip301_streamfunction) s = 1.0_r_def

      lat_pole = pi/2.0_r_def - option(2)*pi
      lon_pole = pi/2.0_r_def + option(3)*pi
 
      psi(1) = 0.0_r_def
      psi(2) = 0.0_r_def
      psi(3) = s * option(1) * scaled_radius * &
               (sin(lat_pole)*sin(chi(2)) - cos(lat_pole)*sin(chi(1)-lon_pole)*cos(chi(2)) )
    case default
      write( log_scratch_space, '(A)' )  'Invalid streamfunction profile choice, stopping'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

end function analytic_streamfunction

end module analytic_streamfunction_profiles_mod
