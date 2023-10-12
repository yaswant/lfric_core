!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Kernel to do a vertical semi-Lagragian advection of rho-type field
!!        using the vertical wind (w) only.
!> @Details The 1D vertical advective transport equation for a W3 variable
!!          is solved using a semi-Lagragian advection scheme.

module vertical_sl_rho_kernel_mod

use argument_mod,          only : arg_type,              &
                                  GH_FIELD, GH_SCALAR,   &
                                  GH_REAL, GH_INTEGER,   &
                                  GH_READWRITE, GH_READ, &
                                  CELL_COLUMN, GH_LOGICAL
use fs_continuity_mod,     only : W2v, W3, Wtheta
use constants_mod,         only : r_tran, i_def, l_def, EPS_R_TRAN
use kernel_mod,            only : kernel_type
! TODO #3011: these config options should be passed through as arguments
use transport_config_mod,  only : vertical_sl_order_cubic,   &
                                  vertical_sl_order_quintic, &
                                  vertical_sl_order_cubic_hermite
use transport_enumerated_types_mod, only : vertical_monotone_none,           &
                                           vertical_monotone_strict,         &
                                           vertical_monotone_relaxed,        &
                                           vertical_monotone_order_constant, &
                                           vertical_monotone_order_linear,   &
                                           vertical_monotone_order_high
use vertical_sl_theta_kernel_mod,   only : monotone_cubic_sl,      &
                                           monotone_quintic_sl,    &
                                           compute_cubic_coeffs,   &
                                           compute_quintic_coeffs, &
                                           compute_cubic_hermite_coeffs
implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed
!>                                      by the PSy layer.
type, public, extends(kernel_type) :: vertical_sl_rho_kernel_type
  private
  type(arg_type) :: meta_args(7) = (/                      &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W2v), & ! departure points
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),  & ! rho
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),  & ! theta-height
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),           & ! sl-order
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),           & ! monotone scheme
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),           & ! monotone order
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)            & ! log_space
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: vertical_sl_rho_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: vertical_sl_rho_code

contains

!-------------------------------------------------------------------------------
!> @details This kernel calculates the departure point of w/theta-points using
!!          only w (i.e., vertical motion only), then interpolate theta at the
!!          departure point using 1d-Cubic-Lagrange interpolation.
!> @param[in]     nlayers      The number of layers
!> @param[in]     dep_pts_z    The vertical departure distance used for SL advection
!> @param[in,out] rho          The rho field at time level n --> rho after SL-advection
!> @param[in]     theta_height            The height of theta-points
!> @param[in]     sl_order                Order of the SL scheme
!> @param[in]     vertical_monotone       The monotone scheme
!> @param[in]     vertical_monotone_order Order of the monotone scheme
!> @param[in]     log_space    Switch to use natural logarithmic space
!!                             for the SL interpolation
!> @param[in]     ndf_w2v      The number of degrees of freedom per cell
!!                             on W2v space
!> @param[in]     undf_w2v     The number of unique degrees of freedom
!!                             on W2v space
!> @param[in]     map_w2v      The dofmap for the cell at the base of the column
!!                             on W2v space
!> @param[in]     ndf_w3       The number of degrees of freedom per cell
!!                             on w3 space
!> @param[in]     undf_w3      The number of unique degrees of freedom
!!                             on w3 space
!> @param[in]     map_w3       The dofmap for the cell at the base of the column
!!                             on w3 space
!-------------------------------------------------------------------------------

subroutine vertical_sl_rho_code( nlayers,                            &
                                 dep_pts_z,                          &
                                 rho,                                &
                                 theta_height,                       &
                                 sl_order,                           &
                                 vertical_monotone,                  &
                                 vertical_monotone_order,            &
                                 log_space,                          &
                                 ndf_w2v, undf_w2v, map_w2v,         &
                                 ndf_w3, undf_w3, map_w3,            &
                                 ndf_wtheta, undf_wtheta, map_wtheta )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                         :: nlayers
  integer(kind=i_def), intent(in)                         :: ndf_w2v
  integer(kind=i_def), intent(in)                         :: undf_w2v
  integer(kind=i_def), dimension(ndf_w2v), intent(in)     :: map_w2v
  integer(kind=i_def), intent(in)                         :: ndf_w3
  integer(kind=i_def), intent(in)                         :: undf_w3
  integer(kind=i_def), dimension(ndf_w3), intent(in)      :: map_w3
  integer(kind=i_def), intent(in)                         :: ndf_wtheta
  integer(kind=i_def), intent(in)                         :: undf_wtheta
  integer(kind=i_def), dimension(ndf_wtheta), intent(in)  :: map_wtheta
  real(kind=r_tran), dimension(undf_w2v), intent(in)      :: dep_pts_z
  real(kind=r_tran), dimension(undf_w3), intent(inout)    :: rho
  real(kind=r_tran), dimension(undf_wtheta),   intent(in) :: theta_height
  integer(kind=i_def), intent(in)  :: sl_order, vertical_monotone,  &
                                      vertical_monotone_order
  logical(kind=l_def), intent(in)  :: log_space
  !
  ! locals
  !
  real(kind=r_tran),  dimension(nlayers)    :: zm, zmd, f0, fd
  real(kind=r_tran),  dimension(nlayers+1)  :: dist, zl, zld
  integer(kind=i_def)                       :: k, nz, nzl, km1, km2, si
  real(kind=r_tran)                         :: d, r, sr
  real(kind=r_tran),   dimension(nlayers,4) :: cc
  integer(kind=i_def), dimension(nlayers,4) :: sc
  real(kind=r_tran),   dimension(nlayers,2) :: cl
  real(kind=r_tran),   dimension(nlayers)   :: dz
  real(kind=r_tran),   dimension(nlayers,6) :: cq
  integer(kind=i_def), dimension(nlayers,6) :: sq

  nz  = nlayers
  nzl = nz + 1_i_def

  ! Extract and fill local column from global data
  ! Map departure points into 1d-array dist
  ! Map theta-height to 1d-zl (zl is cell-edges in the vertical)
  !
  do k = 0, nlayers
    dist(k+1) = dep_pts_z(map_w2v(1)+k)
      zl(k+1) = theta_height(map_wtheta(1)+k)
  end do
  !Map global field into 1d-array f0
  do k = 0, nlayers - 1
    f0(k+1) = rho(map_w3(1)+k)
  end do

  ! Apply log to f0 if required
  ! If using the log_space option, f0 is forced to be positive
  if (log_space) then
    do k = 1, nlayers
      f0(k) = log(max(EPS_R_TRAN,abs(f0(k))))
    end do
  end if

  ! Recover the physical departure points of cell edges zld
  do k = 1, nzl
     d     = abs(dist(k))
     sr    = sign(1.0_r_tran,dist(k))
     si    = int(sr,i_def)
     km1   = int( d,i_def)
     r     = d - real(km1,r_tran)
     km1   = k - km1*si
     km2   = km1 - si
     km1   = max(1_i_def, min(km1,nzl))
     km2   = max(1_i_def, min(km2,nzl))
     zld(k) = zl(km1) - sr*r*abs(zl(km2)-zl(km1))
     zld(k) = min(zl(nzl),max(zl(1),zld(k)))
  end do

  ! Create a local 1d SL problem for cell centres
  ! Compute zm from zl and zmd from zld of cells where,
  ! (zm,zmd) are averaged from (zl,zld)
  !
  do k = 1, nz
    zm(k) = 0.5_r_tran*(zl(k) +  zl(k+1))
  end do
  do k = 1, nz
    zmd(k) = 0.5_r_tran*( zld(k) + zld(k+1) )
    zmd(k) = min(zm(nz),max(zm(1),zmd(k)))
  end do
  !Define the spacing dz between zm-points
  do k = 1, nz - 1
    dz(k) = zm(k+1) - zm(k)
  end do
  dz(nz) = dz(nz - 1)

  select case( sl_order )
  case ( vertical_sl_order_cubic_hermite )

    ! Compute the cubic-hermite interpolation coefficients

    call compute_cubic_hermite_coeffs(zmd,zm,dz,sc,cc,cl,nz,nz)

    ! Do field-interpolation

    do k = 1, nz
      fd(k) = cc(k,1)*f0(sc(k,1)) + cc(k,2)*f0(sc(k,2)) + &
              cc(k,3)*f0(sc(k,3)) + cc(k,4)*f0(sc(k,4))
    end do

    ! If using log_space then convert back
    if (log_space) then
      do k = 1, nz
        fd(k) = exp(fd(k))
        f0(k) = exp(f0(k))
      end do
    end if

    !
    ! Enforce monotonicity if required
    !
    if ( vertical_monotone /= vertical_monotone_none ) then
      ! Apply monotonicity
      call monotone_cubic_sl(fd,f0,sc,cl,vertical_monotone, &
                              vertical_monotone_order,1,nz  )
    end if

  case ( vertical_sl_order_cubic )

    ! Compute the cubic-interpolation coefficients

    call compute_cubic_coeffs(zmd,zm,dz,sc,cc,cl,nz,nz)

    ! Do field-interpolation

    do k = 1, nz
      fd(k) = cc(k,1)*f0(sc(k,1)) + cc(k,2)*f0(sc(k,2)) + &
              cc(k,3)*f0(sc(k,3)) + cc(k,4)*f0(sc(k,4))
    end do

    ! If using log_space then convert back
    if (log_space) then
      do k = 1, nz
        fd(k) = exp(fd(k))
        f0(k) = exp(f0(k))
      end do
    end if

    !
    ! Enforce monotonicity if required
    !
    if ( vertical_monotone /= vertical_monotone_none ) then
      ! Apply monotonicity
      call monotone_cubic_sl(fd,f0,sc,cl,vertical_monotone, &
                              vertical_monotone_order,1,nz  )
    end if

  case ( vertical_sl_order_quintic )

    ! Compute the quintic-interpolation coefficients

    call compute_quintic_coeffs(zmd,zm,dz,sq,cq,cl,nz,nz)

    ! Do field-interpolation

    do k = 1, nz
      fd(k) = cq(k,1)*f0(sq(k,1)) + cq(k,2)*f0(sq(k,2)) + &
              cq(k,3)*f0(sq(k,3)) + cq(k,4)*f0(sq(k,4)) + &
              cq(k,5)*f0(sq(k,5)) + cq(k,6)*f0(sq(k,6))
    end do

    ! If using log_space then convert back
    if (log_space) then
      do k = 1, nz
        fd(k) = exp(fd(k))
        f0(k) = exp(f0(k))
      end do
    end if

    ! Enforce monotonicity if required

    if ( vertical_monotone /= vertical_monotone_none ) then
      ! Apply monotonicity
      call monotone_quintic_sl(fd,f0,sq,cl,vertical_monotone, &
                               vertical_monotone_order,1,nz   )
    end if
  end select

  ! Remap the column answer back to the global data

  do k=0,nlayers - 1
    rho(map_w3(1)+k) = fd(k+1)
  end do

end subroutine vertical_sl_rho_code

end module vertical_sl_rho_kernel_mod
