!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel computes the geopotential field

!> @detail Computes the geopotential field Phi = g*r or g*z for Cartesian
!!         domains

module compute_geopotential_kernel_mod

use argument_mod,         only : arg_type, func_type,                     &
                                 GH_FIELD, GH_READ, GH_WRITE,             &
                                 W0, GH_BASIS,                            &
                                 CELLS
use base_mesh_config_mod, only : geometry, &
                                 base_mesh_geometry_spherical
use constants_mod,        only : r_def
use coord_transform_mod,  only : xyz2llr
use kernel_mod,           only : kernel_type
use planet_config_mod,    only : gravity

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: compute_geopotential_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_FIELD,   GH_WRITE, W0),                             &
       arg_type(GH_FIELD*3, GH_READ, W0)                               &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass :: compute_geopotential_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface compute_geopotential_kernel_type
   module procedure compute_geopotential_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public compute_geopotential_code
contains

type(compute_geopotential_kernel_type) function compute_geopotential_kernel_constructor() result(self)
  return
end function compute_geopotential_kernel_constructor

!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf The number of degrees of freedom per cell
!! @param[in] undf The number of unique degrees of freedom
!! @param[in] map Integer array holding the dofmap for the cell at the base of the column
!! @param[inout] theta Real array, the actual data
!! @param[in] chi_1 Real array, the physical x coordinates
!! @param[in] chi_2 Real array, the physical y coordinates
!! @param[in] chi_3 Real array, the physical z coordinates
subroutine compute_geopotential_code(nlayers,phi, &
                                     chi_1,chi_2,chi_3, &
                                     ndf,undf,map)

  !Arguments
  integer, intent(in) :: nlayers, ndf, undf
  integer, dimension(ndf), intent(in) :: map
  real(kind=r_def), dimension(undf), intent(inout)   :: phi
  real(kind=r_def), dimension(undf), intent(in)    :: chi_1
  real(kind=r_def), dimension(undf), intent(in)    :: chi_2
  real(kind=r_def), dimension(undf), intent(in)    :: chi_3

  !Internal variables
  integer          :: df, k
  real(kind=r_def) :: x(3)
  real(kind=r_def) :: lat, lon, r

  if ( geometry == base_mesh_geometry_spherical ) then
    do k = 0, nlayers-1
      do df = 1, ndf
         x(1) = chi_1(map(df) + k)
         x(2) = chi_2(map(df) + k)
         x(3) = chi_3(map(df) + k)
         call xyz2llr(x(1), x(2), x(3), lon, lat, r)

         phi(map(df) + k) =  gravity*r
       end do
    end do
  else
    do k = 0, nlayers-1
      do df = 1, ndf
         phi(map(df) + k) = gravity*chi_3(map(df) + k)
      end do
    end do
  end if

end subroutine compute_geopotential_code

end module compute_geopotential_kernel_mod
