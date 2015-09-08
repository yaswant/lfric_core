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
use kernel_mod,              only : kernel_type
use argument_mod,            only: arg_type,                  &
                                   GH_FIELD, GH_WRITE, GH_READ, &
                                   W0,                        &
                                   CELLS
use constants_mod,           only: PI, r_def, earth_radius, L_NONLINEAR, &
                                   L_COLD_BUBBLE, L_GRAVITY_WAVE
use coord_transform_mod,     only: xyz2llr
use slush_mod,               only: l_spherical
use generate_global_gw_fields_mod, only: generate_global_gw_pert
use reference_profile_mod,   only: reference_profile
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
  ! compute the pointwise theta profile
  if ( L_NONLINEAR ) then 
    nl = 1.0
  else
    nl = 0.0
  end if
  if ( l_spherical ) then
    do k = 0, nlayers-1
      do df = 1, ndf
         x(1) = chi_1(map(df) + k)
         x(2) = chi_2(map(df) + k)
         x(3) = chi_3(map(df) + k)
         call reference_profile(exner_ref, rho_ref, theta_ref, x)         
         call xyz2llr(x(1), x(2), x(3), lon, lat, r)
         theta_pert = generate_global_gw_pert(lon,lat,r-earth_radius)

         theta(map(df) + k) =  theta_pert + nl*theta_ref
       end do
    end do
  else
    do k = 0, nlayers-1
      do df = 1, ndf
         x(1) = chi_1(map(df) + k)
         x(2) = chi_2(map(df) + k)
         x(3) = chi_3(map(df) + k)
         call reference_profile(exner_ref, rho_ref, theta_ref, x)
        if ( L_GRAVITY_WAVE ) then
          theta(map(df) + k) = THETA0 * sin ( PI * x(3) / H )                        &
                             / ( 1.0_r_def + ( x(1) - XC )**2/A**2 ) + nl*theta_ref
        else
          theta(map(df) + k) = theta_ref

          ! Density current test   
          if ( L_COLD_BUBBLE ) then     
            l = sqrt( ((x(1)-XC)/XR)**2 + ((x(3)-ZC_cold)/ZR)**2 )
            if ( l <= 1.0_r_def ) then
              dt =  15.0_r_def/2.0_r_def*(cos(PI*l)+1.0_r_def)
              theta(map(df) + k) = theta_ref - dt/exner_ref
            end if 
          else
          ! Warm bubble test        
            l = sqrt( ((x(1)-XC))**2 + ((x(3)-ZC_hot))**2 )
            if ( l <= 50.0_r_def ) then
              dt = 0.5_r_def
            else
              dt = 0.5_r_def*exp(-(l-50.0_r_def)**2/(100.0_r_def)**2)
            end if 
            theta(map(df) + k) = theta_ref + dt
          end if
        end if
      end do
    end do
  end if
  
end subroutine initial_theta_code

end module initial_theta_kernel_mod
