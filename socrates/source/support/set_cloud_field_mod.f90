!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Setup cloud fields for Socrates

module set_cloud_field_mod

implicit none
private
public :: set_cloud_field

contains

! @param[in]  nlayers                 Number of layers
! @param[in]  n_profile               Number of atmospheric profiles (columns)
! @param[in]  area_fraction           Large scale cloud area fraction in layers
! @param[in]  cca                     Convective cloud amount (fraction)
! @param[in]  ccw                     Convective in-cloud water (kg/kg)
! @param[in]  t_layer                 Temperature in layers
! @param[in]  latitude                Latitude field
! @param[in]  longitude               Longitude field
! @param[out] i_cloud_representation  Socrates cloud representation option
! @param[out] i_overlap               Socrates cloud overlap option
! @param[out] i_inhom                 Socrates cloud inhomogeneity option
! @param[out] i_drop_re               Socrates droplet effective radius option
! @param[out] rand_seed               Random seed for MCICA cloud generator
! @param[out] n_cloud_layer           Number of cloud layers for Socrates
! @param[out] cloud_frac              Large scale cloud fraction
! @param[out] liq_dim                 Droplet effective radius when supplied
! @param[out] conv_frac               Convective cloud fraction
! @param[out] liq_conv_frac           Convective liquid cloud fraction
! @param[out] ice_conv_frac           Convective ice cloud fraction
! @param[out] liq_conv_mmr            Convective liquid gridbox MMR
! @param[out] ice_conv_mmr            Convective ice gridbox MMR
! @param[out] liq_conv_dim            Convective droplet effective radius
subroutine set_cloud_field(nlayers, n_profile, &
  area_fraction, cca, ccw, t_layer, latitude, longitude, &
  i_cloud_representation, i_overlap, i_inhom, i_drop_re, &
  rand_seed, n_cloud_layer, cloud_frac, liq_dim, conv_frac, &
  liq_conv_frac, ice_conv_frac, liq_conv_mmr, ice_conv_mmr, liq_conv_dim)

use constants_mod, only: r_def, i_def, radians_to_degrees
use conversions_mod, only: zerodegc
use cloud_inputs_mod, only: allicetdegc, starticetkelvin
use radiation_config_mod, only: &
  cloud_representation, cloud_overlap, cloud_inhomogeneity, &
  cloud_representation_no_cloud, &
  cloud_representation_liquid_and_ice, &
  cloud_representation_combined, &
  cloud_representation_conv_strat_liq_ice, &
  cloud_representation_split, &
  cloud_overlap_maximum_random, &
  cloud_overlap_random, &
  cloud_overlap_exponential_random, &
  cloud_inhomogeneity_homogeneous, &
  cloud_inhomogeneity_scaling, &
  cloud_inhomogeneity_mcica, &
  cloud_inhomogeneity_cairns, &
  cloud_inhomogeneity_tripleclouds, &
  droplet_effective_radius, &
  droplet_effective_radius_constant, &
  droplet_effective_radius_liu, &
  constant_droplet_effective_radius
use socrates_runes, only: &
  ip_cloud_representation_off, ip_cloud_representation_ice_water, &
  ip_cloud_representation_combine_ice_water, &
  ip_cloud_representation_csiw, ip_cloud_representation_split_ice_water, &
  ip_overlap_max_random, ip_overlap_random, ip_overlap_exponential_random, &
  ip_inhom_homogeneous, ip_inhom_scaling, ip_inhom_mcica, ip_inhom_cairns, &
  ip_inhom_tripleclouds_2019, &
  ip_droplet_re_default, ip_droplet_re_external, ip_droplet_re_liu
use xios, only: xios_date, xios_get_current_date, &
  xios_date_get_day_of_year, xios_date_get_second_of_day

implicit none

integer(i_def), intent(in) :: nlayers, n_profile
real(r_def),    intent(in) :: area_fraction(nlayers)
real(r_def),    intent(in) :: cca(nlayers), ccw(nlayers)
real(r_def),    intent(in) :: t_layer(nlayers)
real(r_def),    intent(in) :: latitude(n_profile), longitude(n_profile)

integer(i_def), intent(out) :: i_cloud_representation
integer(i_def), intent(out) :: i_overlap, i_inhom
integer(i_def), intent(out) :: i_drop_re
integer(i_def), intent(out) :: rand_seed(n_profile)
integer(i_def), intent(out) :: n_cloud_layer

real(r_def), dimension(nlayers), intent(out) :: cloud_frac, liq_dim
  ! Large-scale cloud fields
real(r_def), dimension(nlayers), intent(out) :: conv_frac, &
  liq_conv_frac, ice_conv_frac, liq_conv_mmr, ice_conv_mmr, liq_conv_dim
  ! Convective cloud fields

type(xios_date) :: datetime
integer(i_def) :: k
real(r_def), parameter :: cloud_frac_min = 0.001_r_def
real(r_def), parameter :: dl = 100.0_r_def
  ! Number of unique seeds per degree of latitude / longitude (for MCICA).
  ! 100 per degree equates to approximately a 1km square area of the
  ! globe (for the Earth).


! Properties of clouds
select case (cloud_representation)
case (cloud_representation_no_cloud)
  i_cloud_representation = ip_cloud_representation_off
case (cloud_representation_liquid_and_ice)
  i_cloud_representation = ip_cloud_representation_ice_water
case (cloud_representation_combined)
  i_cloud_representation = ip_cloud_representation_combine_ice_water
case (cloud_representation_conv_strat_liq_ice)
  i_cloud_representation = ip_cloud_representation_csiw
case (cloud_representation_split)
  i_cloud_representation = ip_cloud_representation_split_ice_water
case default
  i_cloud_representation = ip_cloud_representation_off
end select
select case (cloud_overlap)
case (cloud_overlap_maximum_random)
  i_overlap = ip_overlap_max_random
case (cloud_overlap_random)
  i_overlap = ip_overlap_random
case (cloud_overlap_exponential_random)
  i_overlap = ip_overlap_exponential_random
case default
  i_overlap = ip_overlap_max_random
end select
select case (cloud_inhomogeneity)
case (cloud_inhomogeneity_homogeneous)
  i_inhom = ip_inhom_homogeneous
case (cloud_inhomogeneity_scaling)
  i_inhom = ip_inhom_scaling
case (cloud_inhomogeneity_mcica)
  i_inhom = ip_inhom_mcica
  call xios_get_current_date(datetime)
  ! Generate a unique seed for each gridpoint (actually for each area of
  ! the globe with a size defined by dl) at the given time in seconds.
  rand_seed &
    = int((latitude*radians_to_degrees + 90.0_r_def)*dl, i_def) &
    + int(longitude*radians_to_degrees*dl, i_def)*180_i_def*dl &
    + ((abs(int(datetime%year, i_def)-2000_i_def)*366_i_def) &
    + int(xios_date_get_day_of_year(datetime), i_def))*86400_i_def &
    + int(xios_date_get_second_of_day(datetime), i_def)
case (cloud_inhomogeneity_cairns)
  i_inhom = ip_inhom_cairns
case (cloud_inhomogeneity_tripleclouds)
  i_inhom = ip_inhom_tripleclouds_2019
case default
  i_inhom = ip_inhom_homogeneous
end select
select case (droplet_effective_radius)
case (droplet_effective_radius_constant)
  i_drop_re = ip_droplet_re_external
  do k=1, nlayers
    liq_dim(k) = constant_droplet_effective_radius
    liq_conv_dim(k) = constant_droplet_effective_radius
  end do
case (droplet_effective_radius_liu)
  i_drop_re = ip_droplet_re_liu
case default
  i_drop_re = ip_droplet_re_default
end select

! Treatment of convective cloud
select case (cloud_representation)
case (cloud_representation_combined, &
      cloud_representation_conv_strat_liq_ice, &
      cloud_representation_split)
  do k=1, nlayers
    ! Squeeze large scale cloud into area not filled by convective cloud
    cloud_frac(k) = area_fraction(k) * (1.0_r_def - cca(k))
    ! Convective cloud amount unaltered
    conv_frac(k) = cca(k)
    ! Set liquid fraction of convective cloud using layer temperature
    liq_conv_frac(k) = 1.0_r_def - ( (t_layer(k) - starticeTKelvin) / &
                                     (alliceTdegC-(starticeTKelvin-zerodegc)) )
    liq_conv_frac(k) = MIN(MAX(0.0_r_def, liq_conv_frac(k)), 1.0_r_def)
    ice_conv_frac(k) = 1.0_r_def - liq_conv_frac(k)
    ! Convert liquid and ice fractions to gridbox fractions
    liq_conv_frac(k) = liq_conv_frac(k) * conv_frac(k)
    ice_conv_frac(k) = ice_conv_frac(k) * conv_frac(k)
    ! Calculate conserved gridbox mean liquid and ice MMRs
    liq_conv_mmr(k) = ccw(k) * liq_conv_frac(k)
    ice_conv_mmr(k) = ccw(k) * ice_conv_frac(k)
  end do
case default
  cloud_frac = area_fraction
  conv_frac = 0.0_r_def
end select

! Set the number of cloud layers to the highest cloud layer
n_cloud_layer = 0
do k=nlayers, 1, -1
  if (cloud_frac(k)+conv_frac(k) > cloud_frac_min) then
    n_cloud_layer = k
    exit
  end if
end do

end subroutine set_cloud_field
end module set_cloud_field_mod
