!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Compute the coefficients for reconstructing a
!>        a  1D vertical upwind polynomial representation of a tracer field on the 
!>        faces of a cell
!> @details Compute the coefficients of the flux of a tracer density field using a high order
!>          1D polynomial fit to the integrated tracer values over a given stencil.
!>          The stencil used for the polynomial is centred on the upwind cell for each edge
!>          A symmetric polynomial is used containing all monomials up to the
!>          desired order, i.e. order = 2: 1 + x + x^2
!>          This is exactly fitted over all cells in the stencil
!>          The methodology is inspired by that of Thuburn et.al GMD 2014 for 
!>          2D reconstructions
!>          This method is only valid for lowest order elements 
module poly1d_vert_adv_coeffs_kernel_mod

use argument_mod,      only : arg_type, func_type, mesh_data_type,  &
                              GH_FIELD, GH_INTEGER,                 &
                              GH_WRITE, GH_READ,                    &
                              ANY_SPACE_1,                          &
                              GH_BASIS, CELLS, GH_EVALUATOR
use constants_mod,     only : r_def, i_def, EPS
use fs_continuity_mod, only : Wtheta
use kernel_mod,        only : kernel_type

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: poly1d_vert_adv_coeffs_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                  &
       arg_type(GH_FIELD,   GH_WRITE, Wtheta),                         &
       arg_type(GH_FIELD,   GH_READ,  Wtheta),                         &
       arg_type(GH_FIELD*3, GH_READ,  ANY_SPACE_1),                    &
       arg_type(GH_INTEGER, GH_READ),                                  &
       arg_type(GH_INTEGER, GH_READ)                                   &
       /)
  type(func_type) :: meta_funcs(1) = (/                                &
       func_type(ANY_SPACE_1, GH_BASIS)                                &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: poly1d_vert_adv_coeffs_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! Overload the default structure constructor for function space
interface poly1d_vert_adv_coeffs_kernel_type
   module procedure poly1d_vert_adv_coeffs_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public poly1d_vert_adv_coeffs_code
contains

type(poly1d_vert_adv_coeffs_kernel_type) function poly1d_vert_adv_coeffs_kernel_constructor() result(self)
  implicit none
  return
end function poly1d_vert_adv_coeffs_kernel_constructor

!>@brief Compute the coefficients needed for a 1D vertical reconstruction
!>       of a tracer field on vertical faces
!>@param[in] nlayers Number of vertical layers
!>@param[out] coeff Array of fields to store the coefficients for the polynomial
!!                  reconstruction
!>@param[in]  mdwt Mass matrix diagonal for the W3 space, this is used to give
!!                 the cell volume
!>@param[in] chi1 1st component of the physical coordinate field
!>@param[in] chi2 2nd component of the physical coordinate field
!>@param[in] chi3 3rd component of the physical coordinate field
!>@param[in] ndf_wt Number of degrees of freedom per cell for Wtheta
!>@param[in] undf_wt Total number of degrees of freedom for Wtheta
!>@param[in] map_wt Dofmap of the tracer field
!>@param[in] ndf_wx Number of degrees of freedom per cell for the coordinate space
!>@param[in] undf_wx Total number of degrees of freedom for the coordinate space
!>@param[in] map_wx Dofmap of the coordinate space 
!>@param[in] basis_wx Basis function of the coordinate space evaluated on
!!                    Wtheta nodal points 
!>@param[in] global_order Desired polynomial order for advective computations
!>@param[in] nfaces_v Number of vertical faces (used by PSyclone to size coeff
!!                    array)
subroutine poly1d_vert_adv_coeffs_code(nlayers,                   &
                                       coeff,                     &
                                       mdwt,                      &
                                       chi1, chi2, chi3,          &
                                       ndf_wt,                    &
                                       undf_wt,                   &
                                       map_wt,                    &
                                       ndf_wx,                    &
                                       undf_wx,                   &
                                       map_wx,                    &
                                       basis_wx,                  &
                                       global_order,              &
                                       nfaces_v )


  use matrix_invert_mod,    only: matrix_invert
  use base_mesh_config_mod, only: geometry, &
                                  geometry_spherical
  implicit none
   
  ! Arguments
  integer(kind=i_def), intent(in) :: global_order, nfaces_v
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt, &
                                     ndf_wx, undf_wx

  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), dimension(ndf_wx), intent(in) :: map_wx

  real(kind=r_def), dimension(undf_wt), intent(in) :: mdwt
  real(kind=r_def), dimension(undf_wx), intent(in) :: chi1, chi2, chi3

  real(kind=r_def), dimension(global_order+2, nfaces_v, undf_wt), intent(out) :: coeff

  real(kind=r_def), dimension(1,ndf_wx,ndf_wt), intent(in) :: basis_wx

  ! Local variables
  integer(kind=i_def) :: k, ijk, df, stencil, nmonomial, qp, &
                         m, face, order, kx, spherical, planar
  integer(kind=i_def), allocatable, dimension(:,:)         :: smap
  real(kind=r_def)                                         :: xx, fn, z0, zq, dz
  real(kind=r_def),              dimension(3)              :: x0, xq
  real(kind=r_def), allocatable, dimension(:,:)            :: int_monomial, inv_int_monomial
  real(kind=r_def), allocatable, dimension(:)              :: beta, delta, monomial
  real(kind=r_def),              dimension(global_order+2) :: area

  if ( geometry == geometry_spherical ) then
    spherical = 1.0_r_def
    planar    = 0.0_r_def
  else
    spherical = 0.0_r_def
    planar    = 1.0_r_def
  end if

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

  ! Step 1: Build monomials over all cells in advection stencils

  ! Loop over layers (ignoring first and last points)
  layer_loop: do k = 1, nlayers-1

    ! Compute local order, this is at most global_order but reduces near the 
    ! top and bottom boundary
    order = min(global_order, min(2*(k-1), 2*(nlayers-1-k)))

    ! Number of monomials to use (all polynomials up to total degree of order+1)
    nmonomial = (order + 2)
    allocate( int_monomial(nmonomial, nmonomial),  &
              inv_int_monomial(nmonomial, nmonomial), &
              beta(nmonomial), delta(nmonomial), monomial(nmonomial) )
 
    ! Position vector of bottom of this cell
    x0 = 0.0_r_def
    do df = 1, ndf_wx
      ijk = map_wx( df ) + k
      x0(:) = x0(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*basis_wx(1,df,1)
    end do
    z0 = sqrt(x0(1)**2 + x0(2)**2 + x0(3)**2)*spherical + x0(3)*planar

    ! Position vector of bottom of this cell
    xq = 0.0_r_def
    do df = 1, ndf_wx
      ijk = map_wx( df ) + k
      xq(:) = xq(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*basis_wx(1,df,2)
    end do
    zq = sqrt(xq(1)**2 + xq(2)**2 + xq(3)**2)*spherical + xq(3)*planar

    ! Length of cell to normalise distances by
    dz = zq - z0

    ! Compute the coefficients of each point in the stencil for
    ! each cell when this is the upwind point
    ! Loop over both posivie and negative directions
    face_loop: do face = 0,1
      int_monomial = 0.0_r_def

      ! Loop over all cells in the stencil
      stencil_loop: do stencil = 1, order+2
        area(stencil) = mdwt(map_wt( 1 ) + k + smap(stencil+face,order/2) )
        ! If k + smap(stencil,order/2) == nlayers we need to make sure we use 
        ! coordinate evaluated in cell k = nlayers-1 but with zp=1
        if ( k + smap(stencil+face,order/2) == nlayers ) then
          kx = nlayers-1
          qp = 2
        else
          kx = k + smap(stencil+face,order/2)
          qp = 1
        end if
        ! Evaluate coordinate at theta point of this cell 
        xq = 0.0_r_def
        do df = 1, ndf_wx
          ijk = map_wx( df ) + kx
          xq(:) = xq(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*basis_wx(1,df,qp)
        end do
        zq = sqrt(xq(1)**2 + xq(2)**2 + xq(3)**2)*spherical + xq(3)*planar
 
        ! Second: Compute the local coordinate of each quadrature point from the 
        !         physical coordinate
        xx = (zq - z0)/dz

        ! Third: Compute each needed monomial in terms of the local coordinate
        !        on each quadrature point
        ! Loop over monomials
        do m = 1, nmonomial
          fn = xx**(m-1)
          int_monomial(stencil,m) = int_monomial(stencil,m) + fn*area(stencil)
        end do
      end do stencil_loop

      ! Manipulate the integrals of monomials,
      call matrix_invert(int_monomial,inv_int_monomial,nmonomial)

      ! Initialise polynomial coeffficients to zero
      coeff(:,face+1,map_wt(1)+k) = 0.0_r_def    

      ! Evaluate derivative of polynomial at this theta point
      xx = 0.0_r_def
      monomial(1) = 0.0_r_def
      do stencil = 2, order+2
        monomial(stencil) = real(stencil-1_i_def,r_def)*xx**(stencil-2)
      end do
      do stencil = 1, order+2
        delta(:) = 0.0_r_def
        delta(stencil) = 1.0_r_def
        beta = matmul(inv_int_monomial,delta)
        coeff(stencil,face+1,map_wt(1)+k) = dot_product(monomial,beta)*area(stencil)
      end do
    end do face_loop
    deallocate( int_monomial, inv_int_monomial, beta, delta, monomial )
  end do layer_loop
  ! Enforce zero flux boundary coeffs
  coeff(:,:,map_wt(1) )             = 0.0_r_def
  coeff(:,:,map_wt(2) + nlayers-1 ) = 0.0_r_def

  deallocate( smap )

end subroutine poly1d_vert_adv_coeffs_code

end module poly1d_vert_adv_coeffs_kernel_mod
