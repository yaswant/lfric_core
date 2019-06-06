!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes horizontal tracer values through fitting a high order 1D
!>        upwind reconstruction
!> @details Compute the reconstruction for a tracer field using a high order
!>          polynomial fit to the integrated tracer values. The stencil used for the
!>          polynomial is centred on the upwind cell. 
!>          A 1D cross stencil centred around the cell is used to build the polynomial
!>          representation of the tracer field
!>          This method is only valid for lowest order elements 
module poly1d_adv_recon_kernel_mod

use argument_mod,      only : arg_type, func_type, mesh_data_type,  &
                              GH_FIELD, GH_INTEGER,                 &
                              GH_INC, GH_READ,                      &
                              GH_BASIS, CELLS, GH_EVALUATOR,        &
                              STENCIL, CROSS,                       &
                              reference_element_out_face_normal
use constants_mod,     only : r_def, i_def
use fs_continuity_mod, only : W1, W2, Wtheta
use kernel_mod,        only : kernel_type

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: poly1d_adv_recon_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,   W1),                             &
       arg_type(GH_FIELD,   GH_READ,  W2),                             &
       arg_type(GH_FIELD,   GH_READ,  Wtheta, STENCIL(CROSS)),         &
       arg_type(GH_FIELD,   GH_READ,  Wtheta),                         &
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
  procedure, nopass :: poly1d_adv_recon_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! Overload the default structure constructor for function space
interface poly1d_adv_recon_kernel_type
   module procedure poly1d_adv_recon_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public poly1d_adv_recon_code
contains

  type(poly1d_adv_recon_kernel_type) function poly1d_adv_recon_kernel_constructor() result(self)
  implicit none
  return
end function poly1d_adv_recon_kernel_constructor

!> @brief Computes the horizontal polynomial interpolation of a tracer
!! @param[in]  nlayers Number of layers
!! @param[in,out] reconstruction Recontructed tracer field to compute 
!! @param[in]  wind Wind field
!! @param[in]  tracer Pointwise tracer field to reconstruct
!! @param[in]  coeff Array of polynomial coefficients for interpolation
!! @param[in]  ndf_w1 Number of degrees of freedom per cell
!! @param[in]  undf_w1 Number of unique degrees of freedom for the reconstructed field
!! @param[in]  map_w1 Dofmap for the cell at the base of the column
!! @param[in]  ndf_w2 Number of degrees of freedom per cell
!! @param[in]  undf_w2 Number of unique degrees of freedom for the wind field
!! @param[in]  map_w2 Dofmap for the cell at the base of the column
!! @param[in]  basis_w2 Basis function array evaluated at w1 nodes
!! @param[in]  ndf_wt Number of degrees of freedom per cell
!! @param[in]  undf_wt Number of unique degrees of freedom for the tracer field
!! @param[in]  stencil_size Size of the stencil (number of cells)
!! @param[in]  stencil_map Dofmaps for the stencil
!! @param[in]  order Desired order of polynomial reconstruction
!! @param[in]  nfaces_h Number of horizontal neighbours
!! @param[in]  out_face_normal  Vectors normal to the "out faces" of the
!!                              reference element.
subroutine poly1d_adv_recon_code( nlayers,              &
                                  reconstruction,       &
                                  wind,                 &
                                  tracer,               &
                                  coeff,                &
                                  ndf_w1,               &
                                  undf_w1,              &
                                  map_w1,               &
                                  ndf_w2,               &
                                  undf_w2,              &
                                  map_w2,               &
                                  basis_w2,             &
                                  ndf_wt,               &
                                  undf_wt,              &
                                  stencil_size,         &
                                  stencil_map,          &
                                  order,                &
                                  nfaces_h,             &
                                  out_face_normal )
                                    
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_wt
  integer(kind=i_def), intent(in)                    :: undf_wt
  integer(kind=i_def), intent(in)                    :: ndf_w1
  integer(kind=i_def), intent(in)                    :: undf_w1
  integer(kind=i_def), dimension(ndf_w1), intent(in) :: map_w1
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), intent(in)                    :: order
  integer(kind=i_def), intent(in)                    :: stencil_size
  integer(kind=i_def), intent(in)                    :: nfaces_h

  real(kind=r_def), dimension(undf_w1), intent(inout) :: reconstruction
  real(kind=r_def), dimension(undf_w2), intent(in)    :: wind
  real(kind=r_def), dimension(undf_wt), intent(in)    :: tracer

  real(kind=r_def), dimension(order+1, nfaces_h, undf_wt), intent(in) :: coeff

  real(kind=r_def), dimension(3,ndf_w2,ndf_w1), intent(in) :: basis_w2

  integer(kind=i_def), dimension(ndf_wt,stencil_size), intent(in) :: stencil_map

  real(kind=r_def), intent(in) :: out_face_normal(:,:)

  ! Internal variables
  integer(kind=i_def)                   :: k, df, ij, p, face, stencil
  real(kind=r_def)                      :: direction
  real(kind=r_def), dimension(nfaces_h) :: v_dot_n
  real(kind=r_def)                      :: polynomial_tracer

  real(kind=r_def), dimension(order+1,nfaces_h) :: map1d

  ! Compute 1d map from the cross stencil
  ! i.e for order = 2 the stencil map is
  !      | 5 |
  !  | 2 | 1 | 4 |
  !      | 3 |
  ! so map1d is
  ! ( 1, 2, 4 )
  ! ( 1, 3, 5 )
  ! ( 1, 2, 4 )
  ! ( 1, 3, 5 )
  ! First cell is always the centre cell 
  map1d(1,:) = 1
  do face = 1,nfaces_h
    do stencil = 2,order+1
      map1d(stencil,face) = 2 + mod(face+1,2) + 2*(stencil-2)    
    end do
  end do

  do df = 1,nfaces_h
    v_dot_n(df) = dot_product(basis_w2(:,df,df),out_face_normal(:,df))
  end do

  ij = stencil_map(1,1)

  ! Horizontal tracer reconstruction
  ! Bottom point  
  k = 0
  do df = 1,nfaces_h
    ! Check if this is the upwind cell
    direction = wind(map_w2(df) + k )*v_dot_n(df)
    if ( direction >= 0.0_r_def ) then
      polynomial_tracer = 0.0_r_def
      do p = 1,order+1
        stencil = map1d(p,df)
        polynomial_tracer = polynomial_tracer &
                          + tracer( stencil_map(1,stencil) + k )*coeff( p, df, ij+k )
      end do
      reconstruction(map_w1(df) + k ) = polynomial_tracer
    end if
  end do

  do k = 1, nlayers - 1
    do df = 1,nfaces_h
      ! Check if this is the upwind cell
      direction = (wind(map_w2(df) + k ) + wind(map_w2(df) + k-1 ))*v_dot_n(df)
      if ( direction >= 0.0_r_def ) then
        polynomial_tracer = 0.0_r_def
        do p = 1,order+1
          stencil = map1d(p,df)
          polynomial_tracer = polynomial_tracer &
                            + tracer( stencil_map(1,stencil) + k )*coeff( p, df, ij+k )
        end do
        reconstruction(map_w1(df) + k ) = polynomial_tracer
      end if
    end do
  end do
  ! Final point
  k = nlayers - 1
  do df = 1,nfaces_h
    ! Check if this is the upwind cell
    direction = wind(map_w2(df) + k )*v_dot_n(df)
    if ( direction >= 0.0_r_def ) then
      polynomial_tracer = 0.0_r_def
      do p = 1,order+1
        stencil = map1d(p,df)
        polynomial_tracer = polynomial_tracer &
                          + tracer( stencil_map(2,stencil) + k )*coeff( p, df, ij+k+1 )
      end do
      reconstruction(map_w1(df) + k + 1 ) = polynomial_tracer
    end if
  end do

end subroutine poly1d_adv_recon_code

end module poly1d_adv_recon_kernel_mod
