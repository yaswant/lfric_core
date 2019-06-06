!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Compute a stream function psi such that u = curl(psi)

module initial_streamfunc_kernel_mod

use argument_mod,            only : arg_type, func_type,           &
                                    GH_FIELD, GH_INC, GH_READ,     &
                                    ANY_SPACE_9,                   &
                                    GH_BASIS, GH_DIFF_BASIS,       &
                                    CELLS, GH_QUADRATURE_XYoZ
use constants_mod,           only : r_def, i_def, PI
use fs_continuity_mod,       only : W1
use kernel_mod,              only : kernel_type
use initial_wind_config_mod, only : profile, sbr_angle_lat, sbr_angle_lon, u0, v0

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: initial_streamfunc_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  W1),                              &
       ARG_TYPE(GH_FIELD*3, GH_READ, ANY_SPACE_9)                      &
       /)
  type(func_type) :: meta_funcs(2) = (/                                &
       func_type(W1,          GH_BASIS),                               &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                 &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, public, nopass :: initial_streamfunc_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! Overload the default structure constructor for function space
interface initial_streamfunc_kernel_type
   module procedure initial_streamfunc_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public initial_streamfunc_code
contains

!> @todo fparser faile to parse this declaration if the type name has function
!>in it, i.e. initial_streamfunction_kernel_type, when fparser is fixed (issue
!> #198) the kernel names can be changed back
type(initial_streamfunc_kernel_type) function initial_streamfunc_kernel_constructor() result(self)
  implicit none
  return
end function initial_streamfunc_kernel_constructor

!> @brief Computes the righthand side of the galerkin projection of a 
!>        stream function given by an analytical expression
!> @details Computes rhs = int (c . psi dV) for a vector field psi whose
!!          vertical components contain the values of a stream function given by a
!!          chosen analytic expression
!! @param[in] nlayers Number of layers
!! @param[in,out] rhs Right hand side field to compute
!! @param[in] chi_1 X component of the coordinate field
!! @param[in] chi_2 Y component of the coordinate field
!! @param[in] chi_3 Z component of the coordinate field
!! @param[in] ndf Number of degrees of freedom per cell
!! @param[in] undf Total number of degrees of freedom
!! @param[in] map Dofmap for the cell at the base of the column
!! @param[in] basis Basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_chi Number of dofs per cell for the coordinate field
!! @param[in] undf_chi Total number of degrees of freedom
!! @param[in] map_chi Dofmap for the coordinate field
!! @param[in] chi_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] chi_diff_basis Differential of basis functions evaluated at gaussian quadrature points
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine initial_streamfunc_code(nlayers,                            &
                                   rhs,                                &
                                   chi_1, chi_2, chi_3,                &
                                   ndf, undf,                          &
                                   map, basis,                         &
                                   ndf_chi, undf_chi,                  &
                                   map_chi, chi_basis, chi_diff_basis, &
                                   nqp_h, nqp_v, wqp_h, wqp_v          &
                                   )

  use analytic_streamfunction_profiles_mod, only: analytic_streamfunction
  use base_mesh_config_mod,                 only: geometry, &
                                                  geometry_spherical
  use coordinate_jacobian_mod,              only: coordinate_jacobian, &
                                                  coordinate_jacobian_inverse
  use coord_transform_mod,                  only: sphere2cart_vector, xyz2llr

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf, ndf_chi
  integer(kind=i_def), intent(in) :: undf, undf_chi
  integer(kind=i_def), intent(in) :: nqp_h, nqp_v

  integer(kind=i_def), dimension(ndf),     intent(in) :: map
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi

  real(kind=r_def), intent(in), dimension(3,ndf,    nqp_h,nqp_v) :: basis 
  real(kind=r_def), intent(in), dimension(3,ndf_chi,nqp_h,nqp_v) :: chi_diff_basis
  real(kind=r_def), intent(in), dimension(1,ndf_chi,nqp_h,nqp_v) :: chi_basis

  real(kind=r_def), dimension(undf),     intent(inout) :: rhs
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3

  real(kind=r_def), dimension(nqp_h), intent(in) ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) ::  wqp_v

  ! Internal variables
  integer (kind=i_def)                         :: df, k, qp1, qp2
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jacobian, jac_inv
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_cell, chi_2_cell, chi_3_cell
  real(kind=r_def), dimension(3)               :: psi_physical, psi_spherical, xyz, llr
  real(kind=r_def)                             :: integrand
  real(kind=r_def), dimension(2)               :: option2
  real(kind=r_def), dimension(3)               :: option3

  option3 = (/ U0, sbr_angle_lat, sbr_angle_lon /)
  option2 = (/ U0, V0 /)
  do k = 0, nlayers-1
    do df = 1, ndf_chi
      chi_1_cell(df) = chi_1( map_chi(df) + k)
      chi_2_cell(df) = chi_2( map_chi(df) + k)
      chi_3_cell(df) = chi_3( map_chi(df) + k)
    end do
    call coordinate_jacobian(ndf_chi,        &
                             nqp_h,          &
                             nqp_v,          &
                             chi_1_cell,     &
                             chi_2_cell,     &
                             chi_3_cell,     &
                             chi_diff_basis, &
                             jacobian,       &
                             dj)
    call coordinate_jacobian_inverse(nqp_h, nqp_v, jacobian, dj, jac_inv)
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        ! Compute analytical vector streamfunctiontion in physical space
        xyz(:) = 0.0_r_def
        do df = 1, ndf_chi
          xyz(1) = xyz(1) + chi_1_cell(df)*chi_basis(1,df,qp1,qp2)
          xyz(2) = xyz(2) + chi_2_cell(df)*chi_basis(1,df,qp1,qp2)
          xyz(3) = xyz(3) + chi_3_cell(df)*chi_basis(1,df,qp1,qp2)
        end do
        if ( geometry == geometry_spherical ) then
          call xyz2llr(xyz(1), xyz(2), xyz(3), llr(1), llr(2), llr(3))
          psi_spherical = analytic_streamfunction(llr, profile, 3, option3)
          psi_physical = sphere2cart_vector(psi_spherical,llr) 
        else
          psi_physical = analytic_streamfunction(xyz, profile, 2, option2)
        end if
        do df = 1, ndf 
          integrand = dot_product(matmul(transpose(jac_inv(:,:,qp1,qp2)),&
                                         basis(:,df,qp1,qp2)),psi_physical)*dj(qp1,qp2)
          rhs(map(df) + k) = rhs(map(df) + k) &
                           + wqp_h(qp1)*wqp_v(qp2)*integrand
        end do
      end do
    end do
  end do

end subroutine initial_streamfunc_code

end module initial_streamfunc_kernel_mod
