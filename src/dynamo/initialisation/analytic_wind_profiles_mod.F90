!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Functions to compute a variety of analytic profiles
!> @details Collection of functions to return the value of a field at a given
!!          point based upon a specified analytic formula
module analytic_wind_profiles_mod

use constants_mod,      only : r_def, pi
use initial_wind_config_mod, only : &
                               initial_wind_profile_none,                   &
                               initial_wind_profile_solid_body_rotation,    &
                               initial_wind_profile_constant_uv,            &
                               initial_wind_profile_constant_shear_uv,      &
                               initial_wind_profile_dcmip301,               &
                               initial_wind_profile_deep_baroclinic_steady, &
                               initial_wind_profile_deep_baroclinic_perturbed, &
                               initial_wind_profile_vortex

use planet_config_mod,  only : scaled_radius
use log_mod,            only : log_event,                &
                               log_scratch_space,        &
                               LOG_LEVEL_ERROR
use deep_baroclinic_wave_mod, only : deep_baroclinic_wave

implicit none

private :: vortex_wind

public :: analytic_wind

contains

!> @brief Compute a vortex wind field
!> @param[in] lat Latitudinal position
!> @param[in] long Longitudinal position
!> @param[in] radius Radial distance
!> @result u Result wind field vector (u,v,w)
function vortex_wind(lat,long,radius) result(u)
  implicit none
  real(kind=r_def), intent(in)    :: lat
  real(kind=r_def), intent(in)    :: long
  real(kind=r_def), intent(in)    :: radius
  real(kind=r_def), dimension(3)  :: u

  real(kind=r_def) :: lat_pole, lon_pole, lat_dash
  real(kind=r_def) :: r0, v0, V
  real(kind=r_def) :: radial_distance, omega

  ! Equations below have been taken from Nair and Jablonowski, "Moving vortices
  ! on the sphere: A test case for horizontal advection problems", AMS 2008,
  ! equations (20) and (21)
  ! lat_pole and lon_pole denote the position of one of the vortices, the other
  ! vortex lies on the opposite side of the sphere
  ! Parameter values have been taken from the paper and are currently
  ! hard-wired

  lat_pole = 0.0_r_def
  lon_pole = 0.0_r_def

  lat_dash = asin(sin(lat)*sin(lat_pole) + cos(lat)*cos(lat_pole)*cos(long-lon_pole))

  r0 = 3.0_r_def
  v0 = 38.61073731_r_def
  radial_distance = r0*cos(lat_dash)
  V = v0*3.0_r_def*(sqrt(3.0_r_def)/2.0_r_def)*(sinh(radial_distance)/cosh(radial_distance)**3)

  if (abs(radial_distance)<1E-10_r_def) then
    omega = 0.0_r_def
  else
    omega = V/(radius*radial_distance)
  end if

  u(1) = radius*omega*(sin(lat_pole)*cos(lat) - cos(lat_pole)*cos(long-lon_pole)*sin(lat))
  u(2) = radius*omega*cos(lat_pole)*sin(long-lon_pole)
  u(3) = 0.0_r_def
end function vortex_wind

!> @brief Compute an analytic wind field
!> @param[in] chi Position in physical coordinates
!> @param[in] choice Integer defining which specified formula to use
!> @param[in] num_options Number of sclaer options to supply
!> @param[in] option Array of real values used to generate the initial profile
!> @result u Result wind field vector (u,v,w)
function analytic_wind(chi, choice, num_options, option) result(u)

  implicit none

  real(kind=r_def), intent(in) :: chi(3)
  integer,          intent(in) :: choice, num_options
  real(kind=r_def), optional   :: option(num_options)
  real(kind=r_def)             :: u(3)
  real(kind=r_def)             :: s
  real(kind=r_def)             :: pressure, temperature, density

  if ( .not. present(option) ) option(:) = 0.0_r_def

  select case ( choice )

    case ( initial_wind_profile_none )
      u(:) = 0.0_r_def
    case ( initial_wind_profile_solid_body_rotation, & 
           initial_wind_profile_dcmip301)      
      s = 0.5_r_def*(chi(3)/scaled_radius + 1.0_r_def)
      ! Turn off the height variation for the dcmip test
      if ( choice == initial_wind_profile_dcmip301) s = 1.0_r_def 
      u(1) = s * option(1) * ( cos(chi(2))*cos(option(2)*pi) &
                           + sin(chi(1))*sin(chi(2))*sin(option(2)*pi) )
      u(2) = s * option(1) * cos(chi(1))*sin(option(2)*pi)
      u(3) = 0.0_r_def
    case ( initial_wind_profile_constant_uv )
      u(1) = option(1)
      u(2) = option(2)
      u(3) = 0.0_r_def
    case ( initial_wind_profile_constant_shear_uv )
      u(1) = option(1)*chi(3)/option(3)
      u(2) = option(2)*chi(3)/option(3)
      u(3) = 0.0_r_def 
    case ( initial_wind_profile_deep_baroclinic_steady, &
           initial_wind_profile_deep_baroclinic_perturbed)
      call deep_baroclinic_wave(chi(1), chi(2), chi(3)-scaled_radius, &
                                pressure, temperature, density, &
                                u(1), u(2), u(3)) 
    case ( initial_wind_profile_vortex )
      u = vortex_wind(chi(2),chi(1),chi(3))

    case default
      write( log_scratch_space, '(A)' )  'Invalid velocity profile choice, stopping'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

end function analytic_wind  

end module analytic_wind_profiles_mod
