!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel computes the initial theta field for the gravity wave test

!> @detail The kernel computes initial theta perturbation field for the Klemp & Skamarock 
!>         nonhydrostatic gravity wave test

module initial_theta_kernel_mod

use argument_mod,                  only: arg_type,                          &
                                         GH_FIELD, GH_WRITE, GH_READ,       &
                                         W0,                                &
                                         CELLS
use base_mesh_config_mod,          only: geometry, &
                                         base_mesh_geometry_spherical
use constants_mod,                 only: r_def, PI
use coord_transform_mod,           only: xyz2llr
use formulation_config_mod,        only: nonlinear
use generate_global_gw_fields_mod, only: generate_global_gw_pert
use idealised_config_mod,          only: test,              &
                                         idealised_test_gravity_wave, &
                                         idealised_test_cold_bubble,  &
                                         idealised_test_warm_bubble
use kernel_mod,                    only: kernel_type
use planet_config_mod,             only: scaled_radius
use reference_profile_mod,         only: reference_profile

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: initial_theta_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_FIELD,   GH_WRITE, W0),                             &
       arg_type(GH_FIELD*3, GH_READ,  W0)                              &
       /)
  integer :: iterates_over = CELLS

contains
  procedure, nopass :: initial_theta_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface initial_theta_kernel_type
   module procedure initial_theta_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public initial_theta_code
contains

type(initial_theta_kernel_type) function initial_theta_kernel_constructor() result(self)
  return
end function initial_theta_kernel_constructor

!> @brief The subroutine which is called directly by the psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf The number of degrees of freedom per cell
!! @param[in] undf The number of unique degrees of freedom
!! @param[in] map Integer array holding the dofmap for the cell at the base of the column
!! @param[inout] theta Real array, the actual data
!! @param[in] chi_1 Real array, the physical x coordinates
!! @param[in] chi_2 Real array, the physical y coordinates
!! @param[in] chi_3 Real array, the physical z coordinates
subroutine initial_theta_code(nlayers, &
                              theta, chi_1,chi_2,chi_3, &
                              ndf,undf,map)
  
  !Arguments
  integer, intent(in) :: nlayers, ndf, undf
  integer, dimension(ndf), intent(in) :: map
  real(kind=r_def), dimension(undf), intent(out) :: theta
  real(kind=r_def), dimension(undf), intent(in)    :: chi_1
  real(kind=r_def), dimension(undf), intent(in)    :: chi_2
  real(kind=r_def), dimension(undf), intent(in)    :: chi_3

  !Internal variables
  integer               :: df, k
  real(kind=r_def), parameter :: THETA0 = 0.01_r_def
  real(kind=r_def), parameter :: XC     = 0.0_r_def
  real(kind=r_def), parameter :: YC     = 0.0_r_def
  real(kind=r_def), parameter :: A      = 5000.0_r_def
  real(kind=r_def), parameter :: H      = 10000.0_r_def
  real(kind=r_def)            :: x(3)
  real(kind=r_def)            :: theta_ref, exner_ref, rho_ref
  real(kind=r_def)            :: lat, lon, r
  real(kind=r_def)            :: theta_pert, nl
  real(kind=r_def)            :: l, dt
  real(kind=r_def), parameter :: XR = 4000.0_r_def, &
                                 ZC_cold = 3000.0_r_def, &
                                 ZC_hot = 260.0_r_def, &
                                 ZR = 2000.0_r_def  

  ! Compute the pointwise theta profile
  ! Set up use of nonlinear terms
  if ( nonlinear ) then 
    nl = 1.0
  else
    nl = 0.0
  end if

  ! Set up theta depending on domain shape and choice of idealised test 
  if ( geometry == base_mesh_geometry_spherical ) then   ! SPHERICAL DOMAIN

    ! Gravity wave test only for now
    do k = 0, nlayers-1
      do df = 1, ndf
        x(1) = chi_1(map(df) + k)
        x(2) = chi_2(map(df) + k)
        x(3) = chi_3(map(df) + k)
        ! Calculate reference profile for a chosen idealised test option
        call reference_profile(exner_ref, rho_ref, theta_ref, x, test)
        call xyz2llr(x(1), x(2), x(3), lon, lat, r)
        theta_pert = generate_global_gw_pert(lon,lat,r-scaled_radius)
        theta(map(df) + k) =  theta_pert + nl*theta_ref
      end do
    end do

  else                      ! BIPERIODIC PLANE DOMAIN

    do k = 0, nlayers-1
      do df = 1, ndf
        x(1) = chi_1(map(df) + k)
        x(2) = chi_2(map(df) + k)
        x(3) = chi_3(map(df) + k)
        ! Calculate reference profile for a chosen idealised test option
        call reference_profile(exner_ref, rho_ref, theta_ref, x, test)
        ! Calculate theta a chosen idealised test option
        select case( test )
          case( idealised_test_gravity_wave )  ! Gravity wave test
            theta(map(df) + k) = THETA0 * sin ( PI * x(3) / H )             &
                                 / ( 1.0_r_def + ( x(1) - XC )**2/A**2 ) +  &
                                 nl*theta_ref 
          case( idealised_test_cold_bubble )   ! Density current test
            theta(map(df) + k) = theta_ref     
            l = sqrt( ((x(1)-XC)/XR)**2 + ((x(3)-ZC_cold)/ZR)**2 )
            if ( l <= 1.0_r_def ) then
              dt =  15.0_r_def/2.0_r_def*(cos(PI*l)+1.0_r_def)
              theta(map(df) + k) = theta_ref - dt/exner_ref
            end if 
          case( idealised_test_warm_bubble )   ! Warm bubble test
            theta(map(df) + k) = theta_ref
            l = sqrt( ((x(1)-XC))**2 + ((x(3)-ZC_hot))**2 )
            if ( l <= 50.0_r_def ) then
              dt = 0.5_r_def
            else
              dt = 0.5_r_def*exp(-(l-50.0_r_def)**2/(100.0_r_def)**2)
            end if 
            theta(map(df) + k) = theta_ref + dt
        end select
      end do
    end do

  end if

end subroutine initial_theta_code

end module initial_theta_kernel_mod
