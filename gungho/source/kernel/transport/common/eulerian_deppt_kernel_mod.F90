!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates the Eulerian horizontal departure distance.
!> @details This kernel computes the Eulerian horizontal departure distance
!!          using the advecting departure wind u and time step dt as
!!          \f$\mbox{dist} = x_a - x_d = u dt\f$,
!!          where x_d is the departure point and x_a the arrival point.
!!          This kernel only works with lowest order W2 spaces.

module eulerian_deppt_kernel_mod

  use argument_mod,      only : arg_type,              &
                                GH_FIELD, GH_REAL,     &
                                CELL_COLUMN, GH_WRITE, &
                                GH_READ, GH_SCALAR
  use constants_mod,     only : r_tran, i_def
  use fs_continuity_mod, only : W2, W2h
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: eulerian_deppt_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/               &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W2h), & ! dep pt x
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W2h), & ! dep pt y
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2),  & ! wind
         arg_type(GH_SCALAR, GH_REAL, GH_READ)        & ! dt
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: eulerian_deppt_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: eulerian_deppt_code

contains

  !> @brief Compute the advective increment in x using PPM for the advective fluxes.
  !> @param[in]     nlayers           Number of layers
  !> @param[in,out] dep_pts_x         Departure distance in x
  !> @param[in,out] dep_pts_y         Departure distance in y
  !> @param[in]     wind              Advecting departure wind
  !> @param[in]     dt                Time step
  !> @param[in]     ndf_w2h           Number of degrees of freedom for W2h per cell
  !> @param[in]     undf_w2h          Number of unique degrees of freedom for W2h
  !> @param[in]     map_w2h           Map for W2h
  !> @param[in]     ndf_w2            Number of degrees of freedom for W2 per cell
  !> @param[in]     undf_w2           Number of unique degrees of freedom for W2
  !> @param[in]     map_w2            Map for W2

  subroutine eulerian_deppt_code( nlayers,   &
                                  dep_pts_x, &
                                  dep_pts_y, &
                                  wind,      &
                                  dt,        &
                                  ndf_w2h,   &
                                  undf_w2h,  &
                                  map_w2h,   &
                                  ndf_w2,    &
                                  undf_w2,   &
                                  map_w2 )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w2
    integer(kind=i_def), intent(in) :: ndf_w2
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h

    ! Arguments: Maps
    integer(kind=i_def), dimension(ndf_w2),  intent(in) :: map_w2
    integer(kind=i_def), dimension(ndf_w2h), intent(in) :: map_w2h

    ! Arguments: Fields
    real(kind=r_tran), dimension(undf_w2),  intent(in)    :: wind
    real(kind=r_tran), dimension(undf_w2h), intent(inout) :: dep_pts_x
    real(kind=r_tran), dimension(undf_w2h), intent(inout) :: dep_pts_y
    real(kind=r_tran), intent(in)                         :: dt

    integer(kind=i_def) :: k

    do k = 0,nlayers-1
      ! Eulerian departure distance is just dist = x_a - x_d = u dt
      dep_pts_x( map_w2h(1) + k ) = wind( map_w2(1) + k )*dt
      dep_pts_x( map_w2h(3) + k ) = wind( map_w2(3) + k )*dt
      dep_pts_y( map_w2h(2) + k ) = -wind( map_w2(2) + k )*dt
      dep_pts_y( map_w2h(4) + k ) = -wind( map_w2(4) + k )*dt
    end do

  end subroutine eulerian_deppt_code

end module eulerian_deppt_kernel_mod
