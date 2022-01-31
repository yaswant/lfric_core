!-----------------------------------------------------------------------------
! Copyright (c) 2021,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes the cell integrated mass in the same way as the UM
!>        via a simple rho * dz method.
!>        The value of rho * dz is then multiplied by the square of the
!>        radius from the centre of the Earth to the rho point (a + z_rho)2
!>        and divided  by the square of the radius of the Earth (a2). In this
!>        subroutine, dz will be the distance between theta layers across the
!>        rho level. Only works with lowest order elements.

module compute_column_mass_kernel_mod

  use argument_mod,      only : arg_type, GH_FIELD,          &
                                GH_READWRITE, GH_READ,       &
                                GH_REAL, CELL_COLUMN,        &
                                ANY_DISCONTINUOUS_SPACE_1,   &
                                GH_SCALAR
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W3, WTHETA
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: compute_column_mass_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                         &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),                           &! rho_field
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),                           &! height_w3
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  WTHETA),                       &! height_wth
         arg_type(GH_FIELD,   GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! tot_col_mass
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)                                 &! radius
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: compute_column_mass_code
  end type compute_column_mass_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: compute_column_mass_code

contains

!> @brief Computes the cell integrated mass for a single model column
!! @param[in] nlayers Number of layers
!! @param[in] rho_field Density
!! @param[in,out] tot_col_mass 2D total column mass
!! @param[in] radius Planet radius
!! @param[in] height_w3 Height of density space above surface
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Total number of degrees of freedom for w3
!! @param[in] map_w3 Dof map for the cell at the base of the column for w3
!! @param[in] ndf_wth Number of degrees of freedom per cell for wtheta
!! @param[in] undf_wth Total number of degrees of freedom for wtheta
!! @param[in] map_wth Dof map for the cell at the base of the column for wtheta
!! @param[in] ndf_2d Number of degrees of freedom per cell in 2D field
!! @param[in] undf_2d Total number of degrees of freedom for 2D field
!! @param[in] map_2d Dof map for the cell at the base of the column for 2D field

subroutine compute_column_mass_code(                                   &
                                    nlayers,                           &
                                    rho_field, height_w3,              &
                                    height_wth, tot_col_mass, radius,  &
                                    ndf_w3, undf_w3, map_w3,           &
                                    ndf_wth, undf_wth, map_wth,        &
                                    ndf_2d, undf_2d, map_2d            &
                                   )

  implicit none
  !Arguments
  integer(kind=i_def),                      intent(in)  :: nlayers
  integer(kind=i_def),                      intent(in)  :: ndf_w3, undf_w3
  integer(kind=i_def), dimension(ndf_w3),   intent(in)  :: map_w3
  integer(kind=i_def),                      intent(in)  :: ndf_wth, undf_wth
  integer(kind=i_def), dimension(ndf_wth),  intent(in)  :: map_wth
  integer(kind=i_def),                      intent(in)  :: ndf_2d, undf_2d
  integer(kind=i_def), dimension(ndf_2d),   intent(in)  :: map_2d

  real(kind=r_def),    dimension(undf_w3),  intent(in)    :: rho_field, height_w3
  real(kind=r_def),    dimension(undf_wth), intent(in)    :: height_wth
  real(kind=r_def),    dimension(undf_2d),  intent(inout) :: tot_col_mass
  real(kind=r_def),                         intent(in)    :: radius

  ! Internal variables
  integer(kind=i_def)                          :: k
  ! This is [radius + height_w3]**2 / [radius]**2
  real(kind=r_def)                             :: r2_over_radius2
  real(kind=r_def)                             :: dz_wth
  real(kind=r_def)                             :: mass

  ! Set tot_col_mass to zero
  tot_col_mass( map_2d(1) ) = 0.0_r_def

  do k = 0, nlayers-1

    dz_wth = height_wth( map_wth(1) + k + 1 ) - height_wth( map_wth(1) + k )

    r2_over_radius2 = ((radius + height_w3( map_w3(1) + k ))**2) / radius**2

    mass = ( rho_field( map_w3(1) + k) * dz_wth ) * r2_over_radius2

    tot_col_mass( map_2d(1) ) = tot_col_mass( map_2d(1) ) + mass

  end do

end subroutine compute_column_mass_code

end module compute_column_mass_kernel_mod
