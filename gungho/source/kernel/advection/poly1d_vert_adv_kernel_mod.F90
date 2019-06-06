!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes vertical fluxes through fitting a high order 1D
!>        upwind reconstruction
!> @details Compute the flux for a tracer density field using a high order
!>          polynomial fit to the integrated tracer values. The stencil used for the
!>          polynomial is centred on the upwind cell. 
!>          Near the boundaries the order of reconstruction may be reduced 
!>          if there are not enough points to compute desired order
!>          This method is only valid for lowest order elements 
module poly1d_vert_adv_kernel_mod

use argument_mod,      only : arg_type, func_type, mesh_data_type,  &
                              GH_FIELD, GH_INTEGER,                 &
                              GH_READWRITE, GH_READ,                &
                              GH_BASIS, CELLS, GH_EVALUATOR,        &
                              reference_element_out_face_normal
use constants_mod,     only : r_def, i_def
use fs_continuity_mod, only : W2, Wtheta
use kernel_mod,        only : kernel_type

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: poly1d_vert_adv_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                  &
       arg_type(GH_FIELD,   GH_READWRITE, Wtheta),                     &
       arg_type(GH_FIELD,   GH_READ,      W2),                         &
       arg_type(GH_FIELD,   GH_READ,      Wtheta),                     &
       arg_type(GH_FIELD,   GH_READ,      Wtheta),                     &
       arg_type(GH_INTEGER, GH_READ),                                  &
       arg_type(GH_INTEGER, GH_READ)                                   &
       /)
  type(func_type) :: meta_funcs(1) = (/                                &
       func_type(W2, GH_BASIS)                                         &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass :: poly1d_vert_adv_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! Overload the default structure constructor for function space
interface poly1d_vert_adv_kernel_type
   module procedure poly1d_vert_adv_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public poly1d_vert_adv_code
contains

type(poly1d_vert_adv_kernel_type) function poly1d_vert_adv_kernel_constructor() result(self)
  implicit none
  return
end function poly1d_vert_adv_kernel_constructor

!> @brief Computes the vertical fluxes for a tracer density
!! @param[in]  nlayers Number of layers
!! @param[in,out] advective Advective update to increment 
!! @param[in]  wind Wind field
!! @param[in]  tracer Tracer field to advect
!! @param[in]  coeff Array of polynomial coefficients for interpolation
!! @param[in]  ndf_wt Number of degrees of freedom per cell
!! @param[in]  undf_wt Number of unique degrees of freedom for the tracer field
!! @param[in]  map_wt Cell dofmaps for the tracer space
!! @param[in]  ndf_w2 Number of degrees of freedom per cell
!! @param[in]  undf_w2 Number of unique degrees of freedom for the flux & wind fields
!! @param[in]  map_w2 Dofmap for the cell at the base of the column
!! @param[in]  basis_w2 Basis function array evaluated at wt nodes
!! @param[in]  global_order Desired order of polynomial reconstruction
!! @param[in]  nfaces_v Number of vertical faces (used by PSyclone to size coeff
!!                      array)
subroutine poly1d_vert_adv_code( nlayers,              &
                                 advective,            &
                                 wind,                 &
                                 tracer,               &
                                 coeff,                &
                                 ndf_wt,               &
                                 undf_wt,              &
                                 map_wt,               &
                                 ndf_w2,               &
                                 undf_w2,              &
                                 map_w2,               &
                                 global_order,         &
                                 nfaces_v)
                                    
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_wt
  integer(kind=i_def), intent(in)                    :: undf_wt
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), intent(in)                    :: global_order, nfaces_v

  real(kind=r_def), dimension(undf_wt), intent(inout) :: advective
  real(kind=r_def), dimension(undf_w2), intent(in)    :: wind
  real(kind=r_def), dimension(undf_wt), intent(in)    :: tracer

  real(kind=r_def), dimension(global_order+2, nfaces_v, undf_wt), intent(in) :: coeff

  ! Internal variables
  integer(kind=i_def)            :: k, ij, p, stencil, order, &
                                    m, ijkp, direction
  real(kind=r_def)               :: w 
  real(kind=r_def)               :: polynomial_tracer

  integer(kind=i_def), allocatable, dimension(:,:) :: smap

  ! Compute the offset map for all even orders up to order
  ! The + 3 comes from the minimum number of cells needed
  ! (central cell and the neighbours either side) 
  allocate( smap(global_order+3,0:global_order/2) )
  smap(:,:) = 0
  do m = 0,global_order,2
    do stencil = 1,m+3
      smap(stencil,m/2) = - 1 - m/2 + (stencil-1)
    end do
  end do

  ij = map_wt(1)

  ! Vertical advective update    
  do k = 1, nlayers - 1
    order = min(global_order, min(2*(k-1), 2*(nlayers-1-k)))
    ! Check if this is the upwind cell
    w = sign(1.0_r_def,wind(map_w2(5) + k ))
    direction = 2_i_def - int( 0.5_r_def*(w + abs(w)),i_def )

    polynomial_tracer = 0.0_r_def
    do p = 1,order+2
      ijkp = ij + k + smap(p+direction-1,order/2)
      polynomial_tracer = polynomial_tracer &
                        + tracer( ijkp )*coeff( p, direction, ij+k )
    end do
    advective(map_wt(1) + k ) = advective(map_wt(1) + k ) &
                              + wind(map_w2(5) + k )*polynomial_tracer
      
  end do

  deallocate( smap )
  
end subroutine poly1d_vert_adv_code

end module poly1d_vert_adv_kernel_mod
