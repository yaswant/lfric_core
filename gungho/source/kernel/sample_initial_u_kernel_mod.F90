!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief Kernel computes the rhs for the initialisation of the wind field.

!> @details The kernel computes the rhs of the equation u = u0 where u0 is the
!>          analytically defined wind field. The analytic wind field is projected
!>          onto u using Galerkin projection.

module sample_initial_u_kernel_mod

use argument_mod,            only : arg_type, func_type,           &
                                    GH_FIELD, GH_INC, GH_READ,     &
                                    ANY_SPACE_9, W2,               &
                                    GH_BASIS, GH_DIFF_BASIS,       &
                                    CELLS, GH_QUADRATURE_XYoZ,     &
                                    GH_REAL
use constants_mod,           only : r_def, i_def
use kernel_mod,              only : kernel_type
use initial_wind_config_mod, only : u0, v0

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: sample_initial_u_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  W2),                              &
       ARG_TYPE(GH_FIELD*3, GH_READ, ANY_SPACE_9)                      &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, public, nopass :: sample_initial_u_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! Overload the default structure constructor for function space
interface sample_initial_u_kernel_type
   module procedure sample_initial_u_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

type(sample_initial_u_kernel_type) function sample_initial_u_kernel_constructor() result(self)
  return
end function sample_initial_u_kernel_constructor

!> @brief Compute the right hand side to initialise the wind field.
!! @param[in] nlayers Number of layers
!! @param[in] wind    The velocity vector
!! @param[in] chi_1 X component of the coordinate field
!! @param[in] chi_2 Y component of the coordinate field
!! @param[in] chi_3 Z component of the coordinate field
!! @param[in] ndf Number of degrees of freedom per cell
!! @param[in] undf Total number of degrees of freedom
!! @param[in] map Dofmap for the cell at the base of the column
!! @param[in] ndf_chi Number of dofs per cell for the coordinate field
!! @param[in] undf_chi Total number of degrees of freedom
!! @param[in] map_chi Dofmap for the coordinate field
subroutine sample_initial_u_code(nlayers,                   &
                                 wind,                      & 
                                 chi_1, chi_2, chi_3,       &
                                 ndf, undf, map,            &
                                 ndf_chi, undf_chi, map_chi & 
                                 )

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers, ndf, ndf_chi
  integer(i_def), intent(in) :: undf, undf_chi

  integer(i_def), dimension(ndf),     intent(in) :: map
  integer(i_def), dimension(ndf_chi), intent(in) :: map_chi

  real(r_def), dimension(undf),     intent(inout) :: wind
  real(r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3

  ! Internal variables
  integer(i_def) :: k
  real(r_def)    :: dx, dy, dz, dzdx, dzdy

  do k = 0, nlayers-1
    ! Assumes a linear coordinate field and constant u=U0,v=V0,w=0
    ! u dofs
    dz = 0.5_r_def*(chi_3(map_chi(5)+k) - chi_3(map_chi(1)+k) + chi_3(map_chi(7)+k) - chi_3(map_chi(3)+k))
    dy = 0.5_r_def*(chi_2(map_chi(3)+k) - chi_2(map_chi(1)+k) + chi_2(map_chi(7)+k) - chi_2(map_chi(5)+k))
    wind(map(1)+k) = U0*dy*dz
    dz = 0.5_r_def*(chi_3(map_chi(6)+k) - chi_3(map_chi(2)+k) + chi_3(map_chi(8)+k) - chi_3(map_chi(4)+k))
    dy = 0.5_r_def*(chi_2(map_chi(4)+k) - chi_2(map_chi(2)+k) + chi_2(map_chi(8)+k) - chi_2(map_chi(6)+k))
    wind(map(3)+k) = U0*dy*dz
    ! v dofs
    dz = 0.5_r_def*(chi_3(map_chi(5)+k) - chi_3(map_chi(1)+k) + chi_3(map_chi(6)+k) - chi_3(map_chi(2)+k))
    dx = 0.5_r_def*(chi_1(map_chi(2)+k) - chi_1(map_chi(1)+k) + chi_1(map_chi(6)+k) - chi_1(map_chi(5)+k))
    wind(map(2)+k) = -V0*dx*dz
    dz = 0.5_r_def*(chi_3(map_chi(7)+k) - chi_3(map_chi(3)+k) + chi_3(map_chi(8)+k) - chi_3(map_chi(4)+k))
    dx = 0.5_r_def*(chi_1(map_chi(4)+k) - chi_1(map_chi(3)+k) + chi_1(map_chi(8)+k) - chi_1(map_chi(7)+k))
    wind(map(4)+k) = -V0*dx*dz
    ! W dofs
    dx = 0.5_r_def*(chi_1(map_chi(2)+k) - chi_1(map_chi(1)+k) + chi_1(map_chi(4)+k) - chi_1(map_chi(3)+k))
    dy = 0.5_r_def*(chi_2(map_chi(3)+k) - chi_2(map_chi(1)+k) + chi_2(map_chi(4)+k) - chi_2(map_chi(2)+k))
    dzdx = 0.5_r_def*( chi_3(map_chi(2)+k) - chi_3(map_chi(1)+k) +  chi_3(map_chi(4)+k) - chi_3(map_chi(3)+k) )/dx
    dzdy = 0.5_r_def*( chi_3(map_chi(3)+k) - chi_3(map_chi(1)+k) +  chi_3(map_chi(4)+k) - chi_3(map_chi(2)+k) )/dy
    if ( k > 0)wind(map(5)+k) = -(U0*dzdx+V0*dzdy)*dx*dy
    dx = 0.5_r_def*(chi_1(map_chi(6)+k) - chi_1(map_chi(5)+k) + chi_1(map_chi(8)+k) - chi_1(map_chi(7)+k))
    dy = 0.5_r_def*(chi_2(map_chi(7)+k) - chi_2(map_chi(5)+k) + chi_2(map_chi(8)+k) - chi_2(map_chi(6)+k))
    dzdx = 0.5_r_def*( chi_3(map_chi(6)+k) - chi_3(map_chi(5)+k) +  chi_3(map_chi(8)+k) - chi_3(map_chi(7)+k) )/dx
    dzdy = 0.5_r_def*( chi_3(map_chi(7)+k) - chi_3(map_chi(5)+k) +  chi_3(map_chi(8)+k) - chi_3(map_chi(6)+k) )/dy
    if ( k < nlayers - 1) wind(map(6)+k) = -(U0*dzdx+V0*dzdy)*dx*dy
  end do

end subroutine sample_initial_u_code

end module sample_initial_u_kernel_mod
