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
module poly1d_vert_flux_kernel_mod

use argument_mod,      only : arg_type, func_type, mesh_data_type,  &
                              GH_FIELD, GH_INTEGER,                 &
                              GH_INC, GH_READ,                      &
                              GH_BASIS, CELLS, GH_EVALUATOR,        &
                              reference_element_out_face_normal
use constants_mod,     only : r_def, i_def
use fs_continuity_mod, only : W2, W3
use kernel_mod,        only : kernel_type

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: poly1d_vert_flux_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,   W2),                             &
       arg_type(GH_FIELD,   GH_READ,  W2),                             &
       arg_type(GH_FIELD,   GH_READ,  W3),                             &
       arg_type(GH_FIELD,   GH_READ,  W3),                             &
       arg_type(GH_INTEGER, GH_READ),                                  &
       arg_type(GH_INTEGER, GH_READ)                                   &
       /)
  type(func_type) :: meta_funcs(1) = (/                                &
       func_type(W2, GH_BASIS)                                         &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_EVALUATOR
  type(mesh_data_type) :: meta_init(1) = (/             &
    mesh_data_type( reference_element_out_face_normal ) &
      /)
contains
  procedure, nopass :: poly1d_vert_flux_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! Overload the default structure constructor for function space
interface poly1d_vert_flux_kernel_type
   module procedure poly1d_vert_flux_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public poly1d_vert_flux_code
contains

type(poly1d_vert_flux_kernel_type) function poly1d_vert_flux_kernel_constructor() result(self)
  implicit none
  return
end function poly1d_vert_flux_kernel_constructor

!> @brief Computes the vertical fluxes for a tracer density
!! @param[in]  nlayers Number of layers
!! @param[out] flux Mass flux field to compute 
!! @param[in]  wind Wind field
!! @param[in]  density Tracer density
!! @param[in]  coeff Array of polynomial coefficients for interpolation
!! @param[in]  ndf_w2 Number of degrees of freedom per cell
!! @param[in]  undf_w2 Number of unique degrees of freedom for the flux & wind fields
!! @param[in]  map_w2 Dofmap for the cell at the base of the column
!! @param[in]  basis_w2 Basis function array evaluated at w2 nodes
!! @param[in]  ndf_w3 Number of degrees of freedom per cell
!! @param[in]  undf_w3 Number of unique degrees of freedom for the density field
!! @param[in]  map_w3 Cell dofmaps for the density space
!! @param[in]  global_order Desired order of polynomial reconstruction
!! @param[in]  out_face_normal Vectors normal to the "out faces" of the
!!                             reference element.
!! @param[in]  nfaces_v Number of vertical faces (used by PSyclone to size coeff
!!                     array)
subroutine poly1d_vert_flux_code( nlayers,              &
                                  flux,                 &
                                  wind,                 &
                                  density,              &
                                  coeff,                &
                                  ndf_w2,               &
                                  undf_w2,              &
                                  map_w2,               &
                                  basis_w2,             &
                                  ndf_w3,               &
                                  undf_w3,              &
                                  map_w3,               &
                                  global_order,         &
                                  nfaces_v,             &
                                  out_face_normal )
                                    
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_w3
  integer(kind=i_def), intent(in)                    :: undf_w3
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), intent(in)                    :: global_order, nfaces_v 

  real(kind=r_def), dimension(undf_w2), intent(out)  :: flux
  real(kind=r_def), dimension(undf_w2), intent(in)   :: wind
  real(kind=r_def), dimension(undf_w3), intent(in)   :: density

  real(kind=r_def), dimension(global_order+1, nfaces_v, undf_w3), intent(in) :: coeff

  real(kind=r_def), dimension(3,ndf_w2,ndf_w2), intent(in) :: basis_w2

  real(kind=r_def),    intent(in)  :: out_face_normal(:,:)

  ! Internal variables
  integer(kind=i_def)            :: k, df, ij, p, stencil, order, id, &
                                    m, ijkp
  real(kind=r_def)               :: direction
  real(kind=r_def), dimension(2) :: v_dot_n
  real(kind=r_def)               :: polynomial_density

  integer(kind=i_def), allocatable, dimension(:,:) :: smap

  ! Compute the offset map for all even orders up to order  
  allocate( smap(global_order+1,0:global_order/2) )
  smap(:,:) = 0
  do m = 0,global_order,2
    do stencil = 1,m+1
      smap(stencil,m/2) = - m/2 + (stencil-1)
    end do
  end do

  do df = 5,6    
    v_dot_n(df-4) =  dot_product(basis_w2(:,df,df),out_face_normal(:,df))
  end do
  
  ij = map_w3(1)

  ! Vertical flux computation    
  do k = 0, nlayers - 1
    order = min(global_order, min(2*k, 2*(nlayers-1 - k)))
    do df = 5,6
      id = df - 4 
      ! Check if this is the upwind cell
      direction = wind(map_w2(df) + k )*v_dot_n(id)
      if ( direction >= 0.0_r_def ) then
        polynomial_density = 0.0_r_def
        do p = 1,order+1
          ijkp = ij + k + smap(p,order/2)
          polynomial_density = polynomial_density &
                             + density( ijkp )*coeff( p, id, ij+k )
        end do
        flux(map_w2(df) + k ) = wind(map_w2(df) + k)*polynomial_density
      end if
    end do
  end do

  ! Enforce zero flux boundary conditions
  flux(map_w2(5))           = 0.0_r_def
  flux(map_w2(6)+nlayers-1) = 0.0_r_def

  deallocate( smap )
  
end subroutine poly1d_vert_flux_code

end module poly1d_vert_flux_kernel_mod
