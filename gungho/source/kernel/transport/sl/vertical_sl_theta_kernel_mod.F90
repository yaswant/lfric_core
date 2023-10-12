!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Kernel to do a semi-Lagragian advection of theta-type tracer using
!!        the vertical wind (w) only.
!> @Details The 1D vertical advective transport equation for a Wtheta variable
!!          is solved using a semi-Lagragian advection scheme.

module vertical_sl_theta_kernel_mod

use argument_mod,          only : arg_type,              &
                                  GH_FIELD, GH_SCALAR,   &
                                  GH_READWRITE, GH_READ, &
                                  GH_REAL, GH_INTEGER,   &
                                  CELL_COLUMN, GH_LOGICAL
use fs_continuity_mod,     only : W2v, Wtheta
use constants_mod,         only : r_tran, i_def, l_def, eps, EPS_R_TRAN
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
use vertical_mass_remapping_kernel_mod, only: local_point_1d_array

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed
!>                                      by the PSy layer.
type, public, extends(kernel_type) :: vertical_sl_theta_kernel_type
  private
  type(arg_type) :: meta_args(7) = (/                      &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      W2v   ), & ! departure points
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, Wtheta), & ! theta
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      Wtheta), & ! theta-height
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),           & ! sl-order
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),           & ! monotone scheme
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),           & ! monotone order
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)            & ! log_space
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: vertical_sl_theta_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
 public :: vertical_sl_theta_code
 public :: monotone_cubic_sl
 public :: monotone_quintic_sl
 public :: compute_cubic_coeffs
 public :: compute_quintic_coeffs
 public :: compute_cubic_hermite_coeffs

contains

!-------------------------------------------------------------------------------
!> @details This kernel calculates the departure point of w/theta-points using
!!          only w (i.e., vertical motion only), then interpolates theta at the
!!          departure point using 1d-Cubic-Lagrange interpolation.
!> @param[in]     nlayers      The number of layers
!> @param[in]     dep_pts_z    The vertical departure distance used for SL advection
!> @param[in,out] theta        The theta field at time level n --> theta_d (SL)
!> @param[in]     theta_height            The height of theta-points
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
!> @param[in]     ndf_wtheta   The number of degrees of freedom per cell
!!                             on Wtheta space
!> @param[in]     undf_wtheta  The number of unique degrees of freedom
!!                             on Wtheta space
!> @param[in]     map_wtheta   The dofmap for the cell at the base of the column
!!                             on Wtheta space
!-------------------------------------------------------------------------------

subroutine vertical_sl_theta_code( nlayers,                             &
                                   dep_pts_z,                           &
                                   theta,                               &
                                   theta_height,                        &
                                   sl_order,                            &
                                   vertical_monotone,                   &
                                   vertical_monotone_order,             &
                                   log_space,                           &
                                   ndf_w2v, undf_w2v, map_w2v,          &
                                   ndf_wtheta, undf_wtheta, map_wtheta  )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                        :: nlayers
  integer(kind=i_def), intent(in)                        :: ndf_w2v
  integer(kind=i_def), intent(in)                        :: undf_w2v
  integer(kind=i_def), dimension(ndf_w2v), intent(in)    :: map_w2v
  integer(kind=i_def), intent(in)                        :: ndf_wtheta
  integer(kind=i_def), intent(in)                        :: undf_wtheta
  integer(kind=i_def), dimension(ndf_wtheta), intent(in) :: map_wtheta

  real(kind=r_tran), dimension(undf_w2v), intent(in)       :: dep_pts_z
  real(kind=r_tran), dimension(undf_wtheta), intent(inout) :: theta
  real(kind=r_tran), dimension(undf_wtheta), intent(in)    :: theta_height
  integer(kind=i_def), intent(in)  :: sl_order, vertical_monotone,  &
                                      vertical_monotone_order
  logical(kind=l_def), intent(in)  :: log_space

  real(kind=r_tran)                       :: d, r, sr
  real(kind=r_tran), dimension(nlayers+1) :: zl, zld
  real(kind=r_tran), dimension(nlayers+1) :: dist, theta_local, theta_d_local

  integer(kind=i_def)                         :: k, nzl, km1, km2, si
  real(kind=r_tran),   dimension(nlayers+1,4) :: cc
  integer(kind=i_def), dimension(nlayers+1,4) :: sc
  real(kind=r_tran),   dimension(nlayers+1,2) :: cl
  real(kind=r_tran),   dimension(nlayers+1)   :: dz
  !
  real(kind=r_tran),   dimension(nlayers+1,6) :: cq
  integer(kind=i_def), dimension(nlayers+1,6) :: sq

  ! Extract and fill local column from global variables
  nzl = nlayers + 1_i_def
  do k=0,nlayers
       dist(k+1)        = dep_pts_z(map_w2v(1)+k)
       theta_local(k+1) = theta(map_wtheta(1)+k)
       zl(k+1)          = theta_height(map_wtheta(1)+k)
  end do

  ! Apply log to theta_local if required
  ! If using the log_space option, theta_local is forced to be positive
  if (log_space) then
    do k = 1, nlayers+1
      theta_local(k) = log(max(EPS_R_TRAN,abs(theta_local(k))))
    end do
  end if

  ! Extract physical departure of cell edges
  do k = 1, nzl
    d     = abs(dist(k))
    sr    = sign(1.0_r_tran,dist(k))
    si    = int(sr,i_def)
    km1   = int(d,i_def)
    r     = d - real(km1, r_tran)
    km1   = k - km1*si
    km2   = km1 - si
    km1   = max(1_i_def, min(km1,nzl))
    km2   = max(1_i_def, min(km2,nzl))
    zld(k) = zl(km1) - sr*r*abs(zl(km2)-zl(km1))
    zld(k) = min(zl(nzl),max(zl(1),zld(k)))
  end do
  ! Define the spacing dz between zl-points
  do k = 1, nzl - 1
    dz(k) = zl(k+1) - zl(k)
  end do
  dz(nzl) = dz(nzl-1)

  select case( sl_order )
  case ( vertical_sl_order_cubic_hermite )

    ! Compute the cubic-hermite-interpolation coefficients

    call compute_cubic_hermite_coeffs(zld,zl,dz,sc,cc,cl,nzl,nzl)

    ! Do field-interpolation

    do k = 1, nzl
      theta_d_local(k) = cc(k,1)*theta_local(sc(k,1)) + &
                         cc(k,2)*theta_local(sc(k,2)) + &
                         cc(k,3)*theta_local(sc(k,3)) + &
                         cc(k,4)*theta_local(sc(k,4))
    end do

    ! If using log_space then convert back
    if (log_space) then
      do k = 1, nzl
        theta_d_local(k) = exp(theta_d_local(k))
        theta_local(k)   = exp(theta_local(k))
      end do
    end if

    ! Enforce monotonicity if required

    if ( vertical_monotone /= vertical_monotone_none )  then
      ! Apply monotonicity
      call monotone_cubic_sl(theta_d_local,theta_local,sc,cl,            &
                              vertical_monotone,vertical_monotone_order, &
                              1,nzl                                      )
    end if
  case ( vertical_sl_order_cubic )

    ! Compute the cubic-interpolation coefficients

    call compute_cubic_coeffs(zld,zl,dz,sc,cc,cl,nzl,nzl)

    ! Do field-interpolation

    do k = 1, nzl
      theta_d_local(k) = cc(k,1)*theta_local(sc(k,1)) + &
                         cc(k,2)*theta_local(sc(k,2)) + &
                         cc(k,3)*theta_local(sc(k,3)) + &
                         cc(k,4)*theta_local(sc(k,4))
    end do

    ! If using log_space then convert back
    if (log_space) then
      do k = 1, nzl
        theta_d_local(k) = exp(theta_d_local(k))
        theta_local(k)   = exp(theta_local(k))
      end do
    end if

    ! Enforce monotonicity if required

    if ( vertical_monotone /= vertical_monotone_none )  then
      ! Apply monotonicity
      call monotone_cubic_sl(theta_d_local,theta_local,sc,cl,            &
                              vertical_monotone,vertical_monotone_order, &
                              1,nzl                               )
    end if

  case ( vertical_sl_order_quintic )

   ! Compute the quintic-interpolation coefficients

    call compute_quintic_coeffs(zld,zl,dz,sq,cq,cl,nzl,nzl)

   ! Do field-interpolation

    do k = 1, nzl
      theta_d_local(k) = cq(k,1)*theta_local(sq(k,1)) + &
                         cq(k,2)*theta_local(sq(k,2)) + &
                         cq(k,3)*theta_local(sq(k,3)) + &
                         cq(k,4)*theta_local(sq(k,4)) + &
                         cq(k,5)*theta_local(sq(k,5)) + &
                         cq(k,6)*theta_local(sq(k,6))
    end do

    ! If using log_space then convert back
    if (log_space) then
      do k = 1, nzl
        theta_d_local(k) = exp(theta_d_local(k))
        theta_local(k)   = exp(theta_local(k))
      end do
    end if

    ! Enforce monotonicity if required

    if ( vertical_monotone /= vertical_monotone_none )  then
      ! Apply monotonicity
      call monotone_quintic_sl(theta_d_local,theta_local,sq,cl,    &
                               vertical_monotone,vertical_monotone_order,1,nzl )
    end if
  end select

  ! Remap the column answer back to the global data

  do k=0, nlayers
    theta(map_wtheta(1)+k) = theta_d_local(k+1)
  end do

end subroutine vertical_sl_theta_code
!-------------------------------------------------------------------------------
!> @details This subroutine identify the non-monotone cubic-interpolated values
!!          and modify them using the specified order/option of the modification.
!> @param[inout]  fi        Interpolated field
!> @param[in]     f         The grid-data field
!> @param[in]     sc        The cubic-stencil used for interpolation for each point
!> @param[in]     cl        Linear weights
!> @param[in] vertical_monotone        Option to identify non-monotone points
!> @param[in] vertical_monotone_order  Option of the modification for non-monotone points
!> @param[in]     ns        First index of data-in
!> @param[in]     nf        End index of data-in
!-------------------------------------------------------------------------------
subroutine monotone_cubic_sl(fi,f,sc,cl,vertical_monotone,vertical_monotone_order,ns,nf)
 implicit none
 integer(kind=i_def), intent(in) :: vertical_monotone, vertical_monotone_order, ns, nf
 real(kind=r_tran), dimension(ns:nf), intent(inout)  :: fi(ns:nf)
 real(kind=r_tran), dimension(ns:nf), intent(in)     :: f(ns:nf)
 integer(kind=i_def), dimension(ns:nf,4), intent(in) :: sc
 real(kind=r_tran),   dimension(ns:nf,2), intent(in) :: cl
 ! locals
 integer(kind=i_def), dimension(ns:nf) :: no_mono_id
 integer(kind=i_def)  :: k
 real(kind=r_tran) :: test1,test2,test3,test4,minv,maxv,sigma,xi

 !
 !Identify points/values that are non-monotone using
 !either the strict or relaxed conditions
 !if no_mono_id(:) /= 0, fi is non-monotone
 !   no_mono_id(:) == 1, flat-left  for mod_high
 !   no_mono_id(:) == 2, flat-right for mod_high
 !
 no_mono_id(:) = 0_i_def
 select case( vertical_monotone )
 case ( vertical_monotone_strict )
   do k = ns, nf
     test1 = (fi(k)-f(sc(k,2)))*(f(sc(k,3))-fi(k))
     test2 = (fi(k)-f(sc(k,2)))*(f(sc(k,3))-f(sc(k,2)))
     if ( test1 < 0.0_r_tran ) then
       if (test2  < 0.0_r_tran ) then
         no_mono_id(k) = 1_i_def
       else
         no_mono_id(k) = 2_i_def
       end if
     end if
   end do

 case ( vertical_monotone_relaxed )
   do k = ns, nf
     test1 = (fi(k)-f(sc(k,2)))*(f(sc(k,3))-fi(k))
     test2 = (f(sc(k,2))-f(sc(k,1)))*(f(sc(k,3))-f(sc(k,4)))
     test3 = (fi(k)-f(sc(k,2)))*(f(sc(k,2))-f(sc(k,1)))
     test4 = (fi(k)-f(sc(k,2)))*(f(sc(k,3))-f(sc(k,2)))
     if ( test1 < 0.0_r_tran .and. (test2 < 0.0_r_tran .or. test3 < 0.0_r_tran) ) then
       if (test4  < 0.0_r_tran ) then
         no_mono_id(k) = 1_i_def
       else
         no_mono_id(k) = 2_i_def
       end if
     end if
   end do
 end select
 !
 !Perform modification for identified no-mono points/values
 !
 select case( vertical_monotone_order )
 case ( vertical_monotone_order_constant )
   do k = ns, nf
     if ( no_mono_id(k) /= 0_i_def ) then
       minv = min(f(sc(k,2)), f(sc(k,3)))
       maxv = max(f(sc(k,2)), f(sc(k,3)))
       fi(k) = max( minv, min(fi(k), maxv) )
     end if
   end do

 case ( vertical_monotone_order_linear )
     do k = ns, nf
       if ( no_mono_id(k) /= 0_i_def ) then
         fi(k) = cl(k,1)*f(sc(k,2)) + cl(k,2)*f(sc(k,3))
       end if
     end do

 case ( vertical_monotone_order_high )
     do k = ns, nf
        if (no_mono_id(k) == 1_i_def ) then
          sigma = f(sc(k,3)) -  f(sc(k,2))
          xi = cl(k,2)
          fi(k) =  sigma*xi**3 + f(sc(k,2))
        elseif (no_mono_id(k) == 2_i_def ) then
          sigma = f(sc(k,3)) -  f(sc(k,2))
          xi = cl(k,1)
          fi(k) = -sigma*xi**3 + f(sc(k,3))
        end if
     end do
 end select
end subroutine monotone_cubic_sl

!-------------------------------------------------------------------------------
!> @details This subroutine identify the non-monotone quintic-interpolated values
!!          and modify them using the specified order/option of the modification.
!> @param[inout]  fi        Interpolated field
!> @param[in]     f         The grid-data field
!> @param[in]     sq        The quintic-stencil used for interpolation for each point
!> @param[in]     cl        Linear weights
!> @param[in] vertical_monotone        Option to identify non-monotone points
!> @param[in] vertical_monotone_order  Option of the modification for non-monotone points
!> @param[in]     ns        Start index of data-in
!> @param[in]     nf        End index of data-in
!-------------------------------------------------------------------------------
subroutine monotone_quintic_sl(fi,f,sq,cl,vertical_monotone,vertical_monotone_order,ns,nf)
 implicit none
 integer(kind=i_def), intent(in) :: vertical_monotone,vertical_monotone_order,ns,nf
 real(kind=r_tran), dimension(ns:nf), intent(inout)  :: fi(ns:nf)
 real(kind=r_tran), dimension(ns:nf), intent(in)     :: f(ns:nf)
 integer(kind=i_def), dimension(ns:nf,6), intent(in) :: sq
 real(kind=r_tran),   dimension(ns:nf,2), intent(in) :: cl
 ! locals
 integer(kind=i_def), dimension(ns:nf) :: no_mono_id
 integer(kind=i_def)  :: k
 real(kind=r_tran) :: test1,test2,test3,test4,minv,maxv,sigma, xi

 !
 !Identify points/values that are non-monotone using
 !either the strict or relaxed conditions
 !if no_mono_id(:) # 0, fi is non-monotone
 !   no_mono_id(:) = 1, flat-left  for mod_high
 !   no_mono_id(:) = 2, flat-right for mod_high
 !
 no_mono_id(:) = 0_i_def
 select case( vertical_monotone )
 case ( vertical_monotone_strict )
   do k = ns, nf
      test1 = (fi(k)-f(sq(k,3)))*(f(sq(k,4))-fi(k))
      test2 = (fi(k)-f(sq(k,3)))*(f(sq(k,4))-f(sq(k,3)))
      if ( test1 < 0.0_r_tran ) then
        if (test2  < 0.0_r_tran ) then
          no_mono_id(k) = 1_i_def
        else
          no_mono_id(k) = 2_i_def
        end if
      end if
   end do

 case ( vertical_monotone_relaxed )
   do k = ns, nf
      test1 = (fi(k)-f(sq(k,3)))*(f(sq(k,4))-fi(k))
      test2 = (f(sq(k,3))-f(sq(k,2)))*(f(sq(k,4))-f(sq(k,5)))
      test3 = (fi(k)-f(sq(k,3)) )*(f(sq(k,3))-f(sq(k,2)))
      test4 = (fi(k)-f(sq(k,3)) )*(f(sq(k,4))-f(sq(k,3)))
      if ( test1 < 0.0_r_tran .and. (test2 < 0.0_r_tran .or. test3 < 0.0_r_tran) ) then
         if ( test4  < 0.0_r_tran ) then
           no_mono_id(k) = 1_i_def
         else
           no_mono_id(k) = 2_i_def
         end if
      end if
   end do
 end select
 !
 !Perform modification for identified no-mono points/values
 !
 select case( vertical_monotone_order )
 case ( vertical_monotone_order_constant )
     do k = ns, nf
        if ( no_mono_id(k) /= 0_i_def ) then
           minv = min(f(sq(k,3)), f(sq(k,4)))
           maxv = max(f(sq(k,3)), f(sq(k,4)))
           fi(k) = max( minv, min(fi(k), maxv) )
        end if
     end do

 case ( vertical_monotone_order_linear )
     do k = ns, nf
        if ( no_mono_id(k) /= 0_i_def ) then
           fi(k) = cl(k,1)*f(sq(k,3)) + cl(k,2)*f(sq(k,4))
        end if
     end do

 case ( vertical_monotone_order_high )
     do k = ns, nf
        if ( no_mono_id(k) == 1_i_def ) then
           sigma = f(sq(k,4)) -  f(sq(k,3))
           xi = cl(k,2)
           fi(k) =  sigma*xi**5 + f(sq(k,3))
        elseif ( no_mono_id(k) == 2_i_def ) then
           sigma = f(sq(k,4)) -  f(sq(k,3))
           xi = cl(k,1)
           fi(k) = -sigma*xi**5 + f(sq(k,4))
        end if
     end do
 end select
end subroutine monotone_quintic_sl

!-------------------------------------------------------------------------------
!> @details This subroutine compute cubic-Lagrange weights
!> @param[inout]  zi        The interpolation points
!> @param[in]     zg        The grid points where the data is located
!> @param[in]     sc        The cubic-stencil used for interpolation for each point
!> @param[in]     cc        The cubic interpolation weights
!> @param[in]     cl        Linear weights
!> @param[in]     nzi       The size of interpolation points
!> @param[in]     nzg       The size of the grid-data
!-------------------------------------------------------------------------------
subroutine compute_cubic_coeffs(zi,zg,dz,sc,cc,cl,nzi,nzg)
  implicit none
  integer(kind=i_def),                   intent(in)  :: nzi,nzg
  real(kind=r_tran),   dimension(nzi),   intent(in)  :: zi
  real(kind=r_tran),   dimension(nzg),   intent(in)  :: zg
  real(kind=r_tran),   dimension(nzg),   intent(in)  :: dz
  real(kind=r_tran),   dimension(nzi,4), intent(out) :: cc
  integer(kind=i_def), dimension(nzi,4), intent(out) :: sc
  real(kind=r_tran),   dimension(nzi,2), intent(out) :: cl
  !
  real(kind=r_tran)   :: z1,z2,z3,z4,xi
  real(kind=r_tran)   :: d1,d2,d3,d4
  real(kind=r_tran)   :: n1,n2,n3,n4
  integer(kind=i_def) :: k,km

  do k = 1, nzi
     km = local_point_1d_array(zi(k), zg, 1, nzg)
     xi = (zi(k) - zg(km))/dz(km)

     sc(k,1) = max(1 , km - 1 )
     sc(k,2) = km
     sc(k,3) = min(nzg, km + 1 )
     sc(k,4) = min(nzg, km + 2 )

     z1 = 0.0_r_tran
     z2 = z1 + dz(sc(k,1))
     z3 = z2 + dz(sc(k,2))
     z4 = z3 + dz(sc(k,3))
     xi = z2 + xi*(z3-z2)

     d1 = (z1-z2)*(z1-z3)*(z1-z4)
     d2 = (z2-z1)*(z2-z3)*(z2-z4)
     d3 = (z3-z1)*(z3-z2)*(z3-z4)
     d4 = (z4-z1)*(z4-z2)*(z4-z3)

     n1 = (xi-z2)*(xi-z3)*(xi-z4)
     n2 = (xi-z1)*(xi-z3)*(xi-z4)
     n3 = (xi-z1)*(xi-z2)*(xi-z4)
     n4 = (xi-z1)*(xi-z2)*(xi-z3)
     ! cubic weights
     cc(k,1) = n1/d1
     cc(k,2) = n2/d2
     cc(k,3) = n3/d3
     cc(k,4) = n4/d4
     ! linear weights
     cl(k,1) = (z3-xi)/(z3-z2)
     cl(k,2) = 1.0_r_tran - cl(k,1)

     ! Next to boundaries there is not enough points for cubic
     ! then revert to linear

     if( sc(k,1) == sc(k,2) .or. sc(k,3) == sc(k,4) ) then
        cc(k,1) = 0.0_r_tran
        cc(k,2) = cl(k,1)
        cc(k,3) = cl(k,2)
        cc(k,4) = 0.0_r_tran
     end if
  end do
end subroutine compute_cubic_coeffs

!-------------------------------------------------------------------------------
!> @details This subroutine compute cubic-Hermite weights
!> @param[inout]  zi        The interpolation points
!> @param[in]     zg        The grid points where the data is located
!> @param[in]     sc        The cubic-stencil used for interpolation for each point
!> @param[in]     cc        The cubic interpolation weights
!> @param[in]     cl        Linear weights
!> @param[in]     nzi       The size of interpolation points
!> @param[in]     nzg       The size of the grid-data
!-------------------------------------------------------------------------------
subroutine compute_cubic_hermite_coeffs(zi,zg,dz,sc,cc,cl,nzi,nzg)
  implicit none
  integer(kind=i_def),                   intent(in)  :: nzi,nzg
  real(kind=r_tran),   dimension(nzi),   intent(in)  :: zi
  real(kind=r_tran),   dimension(nzg),   intent(in)  :: zg
  real(kind=r_tran),   dimension(nzg),   intent(in)  :: dz
  real(kind=r_tran),   dimension(nzi,4), intent(out) :: cc
  integer(kind=i_def), dimension(nzi,4), intent(out) :: sc
  real(kind=r_tran),   dimension(nzi,2), intent(out) :: cl
  !
  real(kind=r_tran)    :: xi,alfa,beta,inv_1p_alfa,inv_1p_beta
  real(kind=r_tran)    :: xip2,xip3,c1,c2,c3,c4
  integer(kind=i_def) :: k,km

  do k = 1, nzi
     km = local_point_1d_array(zi(k), zg, 1, nzg)
     xi = (zi(k) - zg(km))/dz(km)

     sc(k,1) = max(1 , km - 1 )
     sc(k,2) = km
     sc(k,3) = min(nzg, km + 1 )
     sc(k,4) = min(nzg, km + 2 )

     alfa = dz(sc(k,1))/dz(sc(k,2))
     beta = dz(sc(k,3))/dz(sc(k,2))
     inv_1p_alfa=1.0_r_tran/(1.0_r_tran + alfa)
     inv_1p_beta=1.0_r_tran/(1.0_r_tran + beta)
     xip2 = xi**2
     xip3 = xi**3
     c1 = 2.0_r_tran*xip3 - 3.0_r_tran*xip2 + 1.0_r_tran
     c2 = 1.0_r_tran - c1
     c3 = xip3 - 2.0_r_tran*xip2 + xi
     c4 = xip3 - xip2
     ! Cubic-hermite weights
     if( sc(k,1) == sc(k,2) ) then
        cc(k,1) = 0.0
        cc(k,2) = c1 - c3 - c4*inv_1p_beta
        cc(k,3) = c2 + c3
        cc(k,4) = c4*inv_1p_beta
     elseif ( sc(k,3) == sc(k,4) ) then
        cc(k,1) = -c3*inv_1p_alfa
        cc(k,2) = c1 - c4
        cc(k,3) = c2 + c4 + c3*inv_1p_alfa
        cc(k,4) = 0.0
     else
        cc(k,1) = -c3*inv_1p_alfa
        cc(k,2) = c1 - c4*inv_1p_beta
        cc(k,3) = c2 + c3*inv_1p_alfa
        cc(k,4) = c4*inv_1p_beta
     end if
     ! linear weights
     cl(k,2) = xi
     cl(k,1) = 1.0_r_tran - cl(k,2)
  end do
end subroutine compute_cubic_hermite_coeffs
!-------------------------------------------------------------------------------
!> @details This subroutine compute quintic-Lagrange weights
!> @param[inout]  zi        The interpolation points
!> @param[in]     zg        The grid points where the data is located
!> @param[in]     sq        The quintic-stencil used for interpolation for each point
!> @param[in]     cq        The quintic interpolation weights
!> @param[in]     cl        Linear weights
!> @param[in]     nzi       The size of interpolation points
!> @param[in]     nzg       The size of the grid-data
!-------------------------------------------------------------------------------
subroutine compute_quintic_coeffs(zi,zg,dz,sq,cq,cl,nzi,nzg)
  implicit none
  integer(kind=i_def),                   intent(in)  :: nzi,nzg
  real(kind=r_tran),   dimension(nzi),   intent(in)  :: zi
  real(kind=r_tran),   dimension(nzg),   intent(in)  :: zg
  real(kind=r_tran),   dimension(nzg),   intent(in)  :: dz
  real(kind=r_tran),   dimension(nzi,6), intent(out) :: cq
  integer(kind=i_def), dimension(nzi,6), intent(out) :: sq
  real(kind=r_tran),   dimension(nzi,2), intent(out) :: cl
  !
  real(kind=r_tran)   :: z1,z2,z3,z4,z5,z6,xi
  real(kind=r_tran)   :: d1,d2,d3,d4,d5,d6
  real(kind=r_tran)   :: n1,n2,n3,n4,n5,n6
  integer(kind=i_def) :: k,km

  do k = 1, nzi
     km = local_point_1d_array(zi(k), zg, 1, nzg)
     xi = (zi(k) - zg(km))/dz(km)

     sq(k,1) = max(1 , km - 2 )
     sq(k,2) = max(1 , km - 1 )
     sq(k,3) = km
     sq(k,4) = min(nzg, km + 1 )
     sq(k,5) = min(nzg, km + 2 )
     sq(k,6) = min(nzg, km + 3 )

     z1 = 0.0_r_tran
     z2 = z1 + dz(sq(k,1))
     z3 = z2 + dz(sq(k,2))
     z4 = z3 + dz(sq(k,3))
     z5 = z4 + dz(sq(k,4))
     z6 = z5 + dz(sq(k,5))
     xi = z3 + xi*(z4-z3)

     d1 = (z1-z2)*(z1-z3)*(z1-z4)*(z1-z5)*(z1-z6)
     d2 = (z2-z1)*(z2-z3)*(z2-z4)*(z2-z5)*(z2-z6)
     d3 = (z3-z1)*(z3-z2)*(z3-z4)*(z3-z5)*(z3-z6)
     d4 = (z4-z1)*(z4-z2)*(z4-z3)*(z4-z5)*(z4-z6)
     d5 = (z5-z1)*(z5-z2)*(z5-z3)*(z5-z4)*(z5-z6)
     d6 = (z6-z1)*(z6-z2)*(z6-z3)*(z6-z4)*(z6-z5)

     if( abs(d1)<eps .or. abs(d2)<eps .or. abs(d3)<eps .or.  &
         abs(d4)<eps .or. abs(d5)<eps .or. abs(d6)<eps ) then
       print*,' moh Q duplicate points z1:z6=',z1,z2,z3,z4,z5,z6
       print*,' moh Q duplicate points z2-z1=',z2-z1
       print*,' moh Q duplicate points z3-z2=',z3-z2
       print*,' moh Q duplicate points z4-z3=',z4-z3
       print*,' moh Q duplicate points z4-z3=',z5-z4
       print*,' moh Q duplicate points z4-z3=',z6-z5
       print*,' moh Q sq(k,:)/nzg = ',sq(k,:),nzg
       print*,' moh Q zg = ',zg
       print*,' moh Q zi/km = ',zi(k),km
     end if

     n1 = (xi-z2)*(xi-z3)*(xi-z4)*(xi-z5)*(xi-z6)
     n2 = (xi-z1)*(xi-z3)*(xi-z4)*(xi-z5)*(xi-z6)
     n3 = (xi-z1)*(xi-z2)*(xi-z4)*(xi-z5)*(xi-z6)
     n4 = (xi-z1)*(xi-z2)*(xi-z3)*(xi-z5)*(xi-z6)
     n5 = (xi-z1)*(xi-z2)*(xi-z3)*(xi-z4)*(xi-z6)
     n6 = (xi-z1)*(xi-z2)*(xi-z3)*(xi-z4)*(xi-z5)

     ! cubic weights
     cq(k,1) = n1/d1
     cq(k,2) = n2/d2
     cq(k,3) = n3/d3
     cq(k,4) = n4/d4
     cq(k,5) = n5/d5
     cq(k,6) = n6/d6
     ! linear weights
     cl(k,1) = (z4-xi)/(z4-z3)
     cl(k,2) = 1.0_r_tran - cl(k,1)
     if( sq(k,1) == sq(k,2) .or. sq(k,5) == sq(k,6) ) then
        ! Revert to cubic weights
          d1 = (z2-z3)*(z2-z4)*(z2-z5)
          d2 = (z3-z2)*(z3-z4)*(z3-z5)
          d3 = (z4-z2)*(z4-z3)*(z4-z5)
          d4 = (z5-z2)*(z5-z3)*(z5-z4)

          n1 = (xi-z3)*(xi-z4)*(xi-z5)
          n2 = (xi-z2)*(xi-z4)*(xi-z5)
          n3 = (xi-z2)*(xi-z3)*(xi-z5)
          n4 = (xi-z2)*(xi-z3)*(xi-z4)

          cq(k,1) = 0.0_r_tran
          cq(k,2) = n1/d1
          cq(k,3) = n2/d2
          cq(k,4) = n3/d3
          cq(k,5) = n4/d4
          cq(k,6) = 0.0_r_tran
      end if
      if( sq(k,2) == sq(k,3) .or. sq(k,4) == sq(k,5) ) then
        ! Revert to linear weights
          cq(k,1) = 0.0_r_tran
          cq(k,2) = 0.0_r_tran
          cq(k,3) = cl(k,1)
          cq(k,4) = cl(k,2)
          cq(k,5) = 0.0_r_tran
          cq(k,6) = 0.0_r_tran
      end if
  end do
end subroutine compute_quintic_coeffs

end module vertical_sl_theta_kernel_mod
