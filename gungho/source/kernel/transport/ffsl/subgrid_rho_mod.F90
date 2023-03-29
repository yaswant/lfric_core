!------------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!------------------------------------------------------------------------------
!> @brief   Routines for calculating coefficients for subgrid rho representation.
!!
!! @details This module contains functions and subroutines which allow quadratic
!!          representation of rho (known as PPM) to be computed.
!------------------------------------------------------------------------------
module subgrid_rho_mod

use constants_mod,                  only: i_def, r_tran, EPS
use transport_enumerated_types_mod, only: horizontal_monotone_none,    &
                                          horizontal_monotone_strict,  &
                                          horizontal_monotone_relaxed, &
                                          vertical_monotone_none,      &
                                          vertical_monotone_strict,    &
                                          vertical_monotone_relaxed

implicit none

private

public :: second_order_vertical_edge
public :: second_order_vertical_gradient
public :: fourth_order_vertical_edge
public :: fourth_order_vertical_edge_strict
public :: fourth_order_vertical_edge_relaxed
public :: horizontal_ppm_coeffs
public :: horizontal_nirvana_coeffs
public :: vertical_nirvana_coeffs
public :: calc_density_at_cell_edge
public :: ppm_output

contains

  !----------------------------------------------------------------------------
  !> @brief Calculates the vertical edge values, taking into account the height
  !!        between layers, using a second-order interpolation.
  !> @details Uses a second-order interpolation to find the vertical cell edge
  !!          values of rho. The vertical grid spacing is used to compute the
  !!          mass, and a quadratic is fit through the cumulative
  !!          mass points. This polynomial is differentiated and evaluated
  !!          at the height of the cell edge, to give the cell edge rho value.
  !!
  !> @param[in]   rho        Density values of two cells which have the ordering
  !!                         | 1 | 2 |
  !> @param[in]   dz         Height of each layer, with index the same as rho
  !> @param[in]   edge_to_do Tells routine which edge to do based on
  !!                         cells       | 1 | 2 |
  !!                         with edges  0   1   2
  !> @param[out]  edge_value The interpolated edge value at edge_to_do
  !----------------------------------------------------------------------------
  subroutine second_order_vertical_edge(rho, dz, edge_to_do, edge_value)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(in)  :: rho(1:2)
    real(kind=r_tran),   intent(in)  :: dz(1:2)
    integer(kind=i_def), intent(in)  :: edge_to_do
    real(kind=r_tran),   intent(out) :: edge_value

    ! Internal Variables
    real(kind=r_tran) :: z(0:2), edge_height
    real(kind=r_tran) :: cmass(0:2)

    ! Get heights of edges starting at 0 for lowest edge in stencil
    z(0) = 0.0_r_tran
    z(1) = z(0) + dz(1)
    z(2) = z(1) + dz(2)

    ! Get edge height to interpolate rho
    edge_height = z(edge_to_do)

    ! Get cumulative mass
    cmass(0) = 0.0_r_tran
    cmass(1) = cmass(0) + dz(1)*rho(1)
    cmass(2) = cmass(1) + dz(2)*rho(2)

    ! Calculate derivative of the quadratic at z = edge_height
    edge_value =   ( 2.0_r_tran*edge_height - z(2) ) / ( z(1) * ( z(1)-z(2) ) ) * cmass(1) &
                 + ( 2.0_r_tran*edge_height - z(1) ) / ( z(2) * ( z(2)-z(1) ) ) * cmass(2)

  end subroutine second_order_vertical_edge

  !----------------------------------------------------------------------------
  !> @brief Calculates the vertical edge gradient, taking into account the height
  !!        between layers, using a second-order method.
  !> @details Uses a second-order method to find the vertical cell edge
  !!          gradient of rho. The vertical grid spacing is used to compute the
  !!          mass, and a quadratic is fit through the cumulative
  !!          mass points. This polynomial is differentiated twice
  !!          to give the cell edge gradient.
  !!
  !> @param[in]   rho           Density values of two cells which have the ordering
  !!                            | 1 | 2 |
  !> @param[in]   dz            Height of each layer, with index the same as rho
  !> @param[out]  edge_gradient The gradient at the edge
  !----------------------------------------------------------------------------
  subroutine second_order_vertical_gradient(rho, dz, edge_gradient)

    implicit none

    ! Arguments
    real(kind=r_tran), intent(in)  :: rho(1:2)
    real(kind=r_tran), intent(in)  :: dz(1:2)
    real(kind=r_tran), intent(out) :: edge_gradient

    ! Internal Variables
    real(kind=r_tran) :: z(0:2)
    real(kind=r_tran) :: cmass(0:2)

    ! Get heights of edges starting at 0 for lowest edge in stencil
    z(0) = 0.0_r_tran
    z(1) = z(0) + dz(1)
    z(2) = z(1) + dz(2)

    ! Get cumulative mass
    cmass(0) = 0.0_r_tran
    cmass(1) = cmass(0) + dz(1)*rho(1)
    cmass(2) = cmass(1) + dz(2)*rho(2)

    ! Calculate second derivative of the quadratic
    edge_gradient =   ( 2.0_r_tran ) / ( z(1) * ( z(1)-z(2) ) ) * cmass(1) &
                    + ( 2.0_r_tran ) / ( z(2) * ( z(2)-z(1) ) ) * cmass(2)

  end subroutine second_order_vertical_gradient

  !----------------------------------------------------------------------------
  !> @brief Calculates the vertical edge values, taking into account the height
  !!        between layers, using a fourth-order interpolation.
  !> @details Uses a fourth-order interpolation to find the vertical cell edge
  !!          values of rho. The vertical grid spacing is used to compute the
  !!          mass, and a high-order polynomial is fit through the cumulative
  !!          mass points. This polynomial is differentiated and evaluated
  !!          at the height of the cell edge, to give the cell edge value.
  !!
  !> @param[in]   rho        Density values of four cells which have the ordering
  !!                         | 1 | 2 | 3 | 4 |
  !> @param[in]   dz         Height of each layer, with index the same as rho
  !> @param[in]   edge_to_do Tells routine which edge to do based on
  !!                         cells       | 1 | 2 | 3 | 4 |
  !!                         with edges  0   1   2   3   4
  !> @param[out]  edge_below The edge value located below layer k
  !!                         (layer k corresponds to cell 3 index above)
  !----------------------------------------------------------------------------
  subroutine fourth_order_vertical_edge(rho, dz, edge_to_do, edge_below)

    implicit none

    real(kind=r_tran),    intent(in)    :: rho(1:4)
    real(kind=r_tran),    intent(in)    :: dz(1:4)
    integer(kind=i_def),  intent(in)    :: edge_to_do
    real(kind=r_tran),    intent(out)   :: edge_below

    real(kind=r_tran) :: z(0:4), dzs(1:4), dzsum, edge_height
    real(kind=r_tran) :: dmass(1:4)
    real(kind=r_tran) :: cmass(0:4)
    real(kind=r_tran) :: poly_mass(1:4)
    real(kind=r_tran) :: dl_dz(1:4)

    integer(kind=i_def) :: i

    ! Get scaling value
    dzsum = sum(dz)

    ! Get scaled dz
    dzs = dz / dzsum

    ! Get heights of edges starting at 0 for lowest edge in stencil
    z(0) = 0.0_r_tran
    do i = 1, 4
      z(i) = z(i-1) + dzs(i)
    end do

    ! Get edge height to interpolate rho to
    edge_height = z(edge_to_do)

    ! Get mass scaled by height
    dmass = rho * dzs

    ! Get cumulative mass
    cmass(0) = 0.0_r_tran
    do i = 1, 4
      cmass(i) = cmass(i-1) + dmass(i)
    end do

    ! Get cumulative mass divided by denominator of polynomial
    poly_mass(1) = cmass(1)/((z(1))*(z(1)-z(2))*(z(1)-z(3))*(z(1)-z(4)))
    poly_mass(2) = cmass(2)/((z(2))*(z(2)-z(1))*(z(2)-z(3))*(z(2)-z(4)))
    poly_mass(3) = cmass(3)/((z(3))*(z(3)-z(1))*(z(3)-z(2))*(z(3)-z(4)))
    poly_mass(4) = cmass(4)/((z(4))*(z(4)-z(1))*(z(4)-z(2))*(z(4)-z(3)))

    ! Calculate derivative of numerator of polynomial at edge height
    dl_dz    = 4.0_r_tran*edge_height**3
    dl_dz(1) = dl_dz(1) - 3.0_r_tran*(z(2)+z(3)+z(4))*edge_height**2 &
               + 2.0_r_tran*(z(3)*z(4) + z(2)*z(3) + z(2)*z(4))*edge_height - z(2)*z(3)*z(4)
    dl_dz(2) = dl_dz(2) - 3.0_r_tran*(z(1)+z(3)+z(4))*edge_height**2 &
               + 2.0_r_tran*(z(3)*z(4) + z(1)*z(3) + z(1)*z(4))*edge_height - z(1)*z(3)*z(4)
    dl_dz(3) = dl_dz(3) - 3.0_r_tran*(z(1)+z(2)+z(4))*edge_height**2 &
               + 2.0_r_tran*(z(2)*z(4) + z(1)*z(2) + z(1)*z(4))*edge_height - z(1)*z(2)*z(4)
    dl_dz(4) = dl_dz(4) - 3.0_r_tran*(z(1)+z(2)+z(3))*edge_height**2 &
               + 2.0_r_tran*(z(2)*z(3) + z(1)*z(2) + z(1)*z(3))*edge_height - z(1)*z(2)*z(3)

    ! Calculate value of edge below layer k
    edge_below = sum( poly_mass * dl_dz )

  end subroutine fourth_order_vertical_edge

  !----------------------------------------------------------------------------
  !> @brief Calculates the vertical edge values, taking into account the height
  !!        between layers, using a fourth-order interpolation, then applying
  !!        strict monotonicity.
  !> @details Uses a fourth-order interpolation to find the vertical cell edge
  !!          values of rho. The vertical grid spacing is used to compute the
  !!          mass, and a high-order polynomial is fit through the cumulative
  !!          mass points. This polynomial is differentiated and evaluated
  !!          at the height of the cell edge, to give the cell edge value.
  !!          Strict monotonicity constraints are then applied.
  !!
  !> @param[in]   rho        Density values of four cells which have the ordering
  !!                         | 1 | 2 | 3 | 4 |
  !> @param[in]   dz         Height of each layer, with index the same as rho
  !> @param[in]   edge_to_do Tells routine which edge to do based on
  !!                         cells       | 1 | 2 | 3 | 4 |
  !!                         with edges  0   1   2   3   4
  !> @param[out]  edge_below The edge value located below layer k
  !!                         (layer k corresponds to cell 3 index above)
  !----------------------------------------------------------------------------
  subroutine fourth_order_vertical_edge_strict(rho, dz, edge_to_do, edge_below)

    implicit none

    real(kind=r_tran),    intent(in)    :: rho(1:4)
    real(kind=r_tran),    intent(in)    :: dz(1:4)
    integer(kind=i_def),  intent(in)    :: edge_to_do
    real(kind=r_tran),    intent(out)   :: edge_below

    real(kind=r_tran) :: t1, tmin, tmax

    ! Get initial unlimited edge value
    call fourth_order_vertical_edge(rho, dz, edge_to_do, edge_below)

    ! Strict Monotonicity
    if ( edge_to_do>0_i_def .AND. edge_to_do<4_i_def) then
      t1 = ( edge_below - rho(edge_to_do) )*( rho(edge_to_do+1) - edge_below )
      if ( t1 < 0.0_r_tran ) then
        tmin = min(rho(edge_to_do+1),rho(edge_to_do))
        tmax = max(rho(edge_to_do+1),rho(edge_to_do))
        edge_below = min( tmax, max(edge_below,tmin) )
      end if
    else if ( edge_to_do == 0_i_def ) then
      edge_below = min( max( rho(2), rho(1) ), max( edge_below, min( rho(2), rho(1) ) ) )
    else if ( edge_to_do == 4_i_def ) then
      edge_below = min( max( rho(4), rho(3) ), max( edge_below, min( rho(4), rho(3) ) ) )
    end if

  end subroutine fourth_order_vertical_edge_strict

  !----------------------------------------------------------------------------
  !> @brief Calculates the vertical edge values, taking into account the height
  !!        between layers, using a fourth-order interpolation, then applying
  !!        relaxed monotonicity.
  !> @details Uses a fourth-order interpolation to find the vertical cell edge
  !!          values of rho. The vertical grid spacing is used to compute the
  !!          mass, and a high-order polynomial is fit through the cumulative
  !!          mass points. This polynomial is differentiated and evaluated
  !!          at the height of the cell edge, to give the cell edge value.
  !!          Relaxed monotonicity constraints are then applied.
  !!
  !> @param[in]   rho        Density values of four cells which have the ordering
  !!                         | 1 | 2 | 3 | 4 |
  !> @param[in]   dz         Height of each layer, with index the same as rho
  !> @param[in]   edge_to_do Tells routine which edge to do based on
  !!                         cells       | 1 | 2 | 3 | 4 |
  !!                         with edges  0   1   2   3   4
  !> @param[out]  edge_below The edge value located below layer k
  !!                         (layer k corresponds to cell 3 index above)
  !----------------------------------------------------------------------------
  subroutine fourth_order_vertical_edge_relaxed(rho, dz, edge_to_do, edge_below)

    implicit none

    real(kind=r_tran),    intent(in)    :: rho(1:4)
    real(kind=r_tran),    intent(in)    :: dz(1:4)
    integer(kind=i_def),  intent(in)    :: edge_to_do
    real(kind=r_tran),    intent(out)   :: edge_below

    real(kind=r_tran) :: t1, t2, t3, tmin, tmax

    ! Get initial unlimited edge value
    call fourth_order_vertical_edge(rho, dz, edge_to_do, edge_below)

    ! Relaxed Monotonicity
    if ( edge_to_do>0_i_def .AND. edge_to_do<4_i_def) then
      t1 = ( edge_below - rho(edge_to_do) )*( rho(edge_to_do+1) - edge_below )
      if ( edge_to_do == 2_i_def ) then
        t2 = ( rho(edge_to_do) - rho(edge_to_do-1) )*( rho(edge_to_do+2) - rho(edge_to_do+1) )
        t3 = ( edge_below - rho(edge_to_do) )*( rho(edge_to_do) - rho(edge_to_do-1) )
      else
        t2 = 1.0_r_tran
        t3 = -1.0_r_tran
      end if
      if ( t1 < 0.0_r_tran .AND. ( t2 >= 0.0_r_tran .OR. t3 <= 0.0_r_tran ) ) then
        tmin = min(rho(edge_to_do+1),rho(edge_to_do))
        tmax = max(rho(edge_to_do+1),rho(edge_to_do))
        edge_below = min( tmax, max(edge_below,tmin) )
      end if
    else if ( edge_to_do == 0_i_def ) then
      edge_below = min( max( rho(2), rho(1) ), max( edge_below, min( rho(2), rho(1) ) ) )
    else if ( edge_to_do == 4_i_def ) then
      edge_below = min( max( rho(4), rho(3) ), max( edge_below, min( rho(4), rho(3) ) ) )
    end if

  end subroutine fourth_order_vertical_edge_relaxed

  !----------------------------------------------------------------------------
  !> @brief  Returns the PPM coefficients (a0, a1, a2) which are a quadratic
  !!         representation of rho within the cell, i.e. rho(x)=a0 + a1*x + a2*x^2
  !!         for 0<=x<=1. The dofmap for the density values is of the form
  !!         | 1 | 2 | 3 | 4 | 5 | where the subgrid coefficients are being
  !!         estimated for cell 3. This is for the horizontal PPM coefficients,
  !!         and assumes uniform grid spacing in the horizontal.
  !!
  !! @param[in]   density        Density values of five cells which have the
  !!                             ordering
  !!                             | 1 | 2 | 3 | 4 | 5 |
  !! @param[out]  coeffs         Coefficients for cell 3 with coeffs(1)=a0,
  !!                             coeffs(2)=a1, coeffs(3)=a2
  !! @param[in]   monotone       Monotone option to ensures no over/undershoots
  !----------------------------------------------------------------------------
  subroutine horizontal_ppm_coeffs(coeffs,density,monotone)

    implicit none

    real(kind=r_tran),    intent(out) :: coeffs(1:3)
    real(kind=r_tran),    intent(in)  :: density(1:5)
    integer(kind=i_def),  intent(in)  :: monotone

    real(kind=r_tran)                :: density_cell_edge_left
    real(kind=r_tran)                :: density_cell_edge_right

    coeffs(:) = 0.0_r_tran

    density_cell_edge_left = calc_density_at_cell_edge(density(1:4),monotone)
    density_cell_edge_right = calc_density_at_cell_edge(density(2:5),monotone)
    call ppm_output(density_cell_edge_left,density_cell_edge_right,density(3),monotone,coeffs)

  end subroutine horizontal_ppm_coeffs

  !----------------------------------------------------------------------------
  !> @brief  Returns the Nirvana coefficients (a0, a1, a2) which are a quadratic
  !!         representation of rho within the cell, i.e. rho(x)=a0 + a1*x + a2*x^2
  !!         for 0<=x<=1. The dofmap for the density values is of the form
  !!         | 1 | 2 | 3 | where the subgrid coefficients are being
  !!         estimated for cell 2. This is for the horizontal Nirvana coefficients,
  !!         and assumes uniform grid spacing in the horizontal. See Leonard et al. (1995).
  !!         The conditions on the coefficients for Nirvana are equivalent to the integral
  !!         of the quadratic equaling the integral of rho in all three cells.
  !!
  !! @param[out]  coeffs    Coefficients for cell 2 with
  !!                        coeffs(1)=a0, coeffs(2)=a1, coeffs(3)=a2
  !! @param[in]   rho       Density values of three cells which have the ordering
  !!                        | 1 | 2 | 3 |
  !! @param[in]   monotone  Monotone option to ensures no over/undershoots
  !----------------------------------------------------------------------------
  subroutine horizontal_nirvana_coeffs(coeffs,rho,monotone)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(out) :: coeffs(1:3)
    real(kind=r_tran),   intent(in)  :: rho(1:3)
    integer(kind=i_def), intent(in)  :: monotone

    ! Internal variables
    real(kind=r_tran) :: t1, t2, t3
    real(kind=r_tran) :: rho_left, rho_right

    ! Initialise coefficients to be zero
    coeffs(:) = 0.0_r_tran

    ! The coefficients are taken from Leonard et al. (1995) for a uniform grid

    coeffs(1) = (-rho(3) + 5.0_r_tran * rho(2) + 2.0_r_tran * rho(1)) / 6.0_r_tran
    coeffs(2) = rho(2) - rho(1)
    coeffs(3) = (rho(3) - 2.0_r_tran * rho(2) + rho(1)) / 2.0_r_tran

    ! Apply monotonicity if needed
    if ( monotone == horizontal_monotone_strict ) then
      ! Strict monotonicity
      t1 = -0.5_r_tran * coeffs(2) / ( coeffs(3) + EPS )
      if ( ( t1 + EPS ) * ( 1.0_r_tran + EPS - t1 ) > 0.0_r_tran ) then
        coeffs(1) = rho(2)
        coeffs(2) = 0.0_r_tran
        coeffs(3) = 0.0_r_tran
      end if
    else if ( monotone == horizontal_monotone_relaxed ) then
      ! Approximate cell edge values
      rho_left  = ( rho(2) + rho(1) ) / 2.0_r_tran
      rho_right = ( rho(2) + rho(3) ) / 2.0_r_tran
      ! Relaxed monotonicity
      t1 = -0.5_r_tran * coeffs(2) / ( coeffs(3) + EPS )
      if ( ( t1 + EPS ) * ( 1.0_r_tran + EPS - t1 ) > 0.0_r_tran ) then
        t2 = ( rho_right - rho(2) ) * ( rho(2) - rho_left )
        t3 = abs( rho(2) - rho_left ) - abs( rho_right - rho(2) )
        if ( t2 < 0.0_r_tran ) then
          coeffs(1) = rho(2)
          coeffs(2) = 0.0_r_tran
          coeffs(3) = 0.0_r_tran
        else
          if ( t3 < 0.0_r_tran ) then
            coeffs(1) = rho_left
            coeffs(2) = 0.0_r_tran
            coeffs(3) = 3.0_r_tran*(rho(2) - rho_left)
          else
            coeffs(1) = -2.0_r_tran*rho_right + 3.0_r_tran*rho(2)
            coeffs(2) =  6.0_r_tran*rho_right - 6.0_r_tran*rho(2)
            coeffs(3) = -3.0_r_tran*rho_right + 3.0_r_tran*rho(2)
          end if
        end if
      end if
    end if

  end subroutine horizontal_nirvana_coeffs

  !----------------------------------------------------------------------------
  !> @brief  Returns the vertical Nirvana coefficients (a0, a1, a2) which are a quadratic
  !!         representation of rho within the cell, i.e. rho(z)=a0 + a1*z + a2*z^2
  !!         for 0<=z<=1. The conditions on the coefficients for Nirvana are equivalent to:
  !!         1) The gradient of the quadratic equaling the gradient of rho at cell edges.
  !!         2) The integral of the quadratic equaling the integral of rho in the cell.
  !!
  !! @param[out]  coeffs      Coefficients for cell 2 with
  !!                          coeffs(1)=a0, coeffs(2)=a1, coeffs(3)=a2
  !! @param[in]   rho         Average density of the cells | 1 | 2 | 3 |
  !! @param[in]   dz          Height of cell 2
  !! @param[out]  grad_below  Estimate of gradient at z = 0 of cell 2, i.e.
  !!                          at the edge between cells 1 and 2
  !! @param[out]  grad_above  Estimate of gradient at z = 1 of cell 2, i.e.
  !!                          at the edge between cells 2 and 3
  !! @param[in]   monotone    Monotone option to ensures no over/undershoots
  !----------------------------------------------------------------------------
  subroutine vertical_nirvana_coeffs(coeffs,rho,dz,grad_below,grad_above,monotone)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(out) :: coeffs(1:3)
    real(kind=r_tran),   intent(in)  :: rho(1:3)
    real(kind=r_tran),   intent(in)  :: dz, grad_below, grad_above
    integer(kind=i_def), intent(in)  :: monotone

    ! Internal Variables
    real(kind=r_tran) :: t1, t2, t3
    real(kind=r_tran) :: rho_left, rho_right

    ! Calculate coefficients
    coeffs(2) = grad_below * dz
    coeffs(3) = ( grad_above * dz - coeffs(2) ) / 2.0_r_tran
    coeffs(1) = rho(2) - coeffs(2) / 2.0_r_tran - coeffs(3) / 3.0_r_tran

    ! Apply monotonicity if needed
    if ( monotone == vertical_monotone_strict ) then
      ! Strict monotonicity
      t1 = -0.5_r_tran * coeffs(2) / ( coeffs(3) + EPS )
      if ( ( t1 + EPS ) * ( 1.0_r_tran + EPS - t1 ) > 0.0_r_tran ) then
        coeffs(1) = rho(2)
        coeffs(2) = 0.0_r_tran
        coeffs(3) = 0.0_r_tran
      end if
    else if ( monotone == vertical_monotone_relaxed ) then
      ! Approximate cell edge values
      rho_left  = ( rho(2) + rho(1) ) / 2.0_r_tran
      rho_right = ( rho(2) + rho(3) ) / 2.0_r_tran
      ! Relaxed monotonicity
      t1 = -0.5_r_tran * coeffs(2) / ( coeffs(3) + EPS )
      if ( ( t1 + EPS ) * ( 1.0_r_tran + EPS - t1 ) > 0.0_r_tran ) then
        t2 = ( rho_right - rho(2) ) * ( rho(2) - rho_left )
        t3 = abs( rho(2) - rho_left ) - abs( rho_right - rho(2) )
        if ( t2 < 0.0_r_tran ) then
          coeffs(1) = rho(2)
          coeffs(2) = 0.0_r_tran
          coeffs(3) = 0.0_r_tran
        else
          if ( t3 < 0.0_r_tran ) then
            coeffs(1) = rho_left
            coeffs(2) = 0.0_r_tran
            coeffs(3) = 3.0_r_tran*(rho(2) - rho_left)
          else
            coeffs(1) = -2.0_r_tran*rho_right + 3.0_r_tran*rho(2)
            coeffs(2) =  6.0_r_tran*rho_right - 6.0_r_tran*rho(2)
            coeffs(3) = -3.0_r_tran*rho_right + 3.0_r_tran*rho(2)
          end if
        end if
      end if
    end if

  end subroutine vertical_nirvana_coeffs

  !----------------------------------------------------------------------------
  !> @brief  Calculates the estimated density at the edge of a cell required for
  !!         using PPM to estimate the quadratic subgrid representation of rho.
  !!         The function is passed four density values from consecutive cells
  !!         (which all lie in the same direction) with the dofmap
  !!         | 1 | 2 | 3 | 4 | and returns the estimated density value between
  !!         cells 2 and 3. The cells are assumed to be uniform in spacing.
  !!         Monotonicity options are provided.
  !!
  !! @param[in]   density            Has dof map of the form | 1 | 2 | 3 | 4 |
  !! @param[in]   monotone           Monotone option to ensures no over/undershoots
  !! @return      density_at_edge    Interpolated density value at edge between
  !!                                 cells 2 and 3.
  !----------------------------------------------------------------------------
  function calc_density_at_cell_edge(density,monotone) result(density_at_edge)

    implicit none

    real(kind=r_tran),   intent(in) :: density(1:4)
    integer(kind=i_def), intent(in) :: monotone

    real(kind=r_tran) :: density_at_edge
    real(kind=r_tran) :: t1, t2, t3, tmax, tmin

    ! As the cell widths are assumed to be constant the edge value reduces to that given in
    ! Colella and Woodward, JCP, 54, 1984, equation (1.9)
    density_at_edge = (7.0_r_tran/12.0_r_tran) * (density(2)+density(3)) &
                     -(1.0_r_tran/12.0_r_tran) * (density(1)+density(4))

    if ( monotone == horizontal_monotone_strict ) then
      ! Strict monotonicity
      t1 = ( density_at_edge - density(2) )*( density(3) - density_at_edge )
      if ( t1 < 0.0_r_tran ) then
         tmin = min(density(3),density(2))
         tmax = max(density(3),density(2))
         density_at_edge = min( tmax, max(density_at_edge,tmin) )
      end if
    else if ( monotone == horizontal_monotone_relaxed ) then
      ! Relaxed monotonicity
      t1 = ( density_at_edge - density(2) )*( density(3) - density_at_edge )
      t2 = ( density(2) - density(1) )*( density(4) - density(3) )
      t3 = ( density_at_edge - density(2) )*( density(2) - density(1) )
      if ( t1 < 0.0_r_tran .AND. ( t2 >= 0.0_r_tran .OR. t3 <= 0.0_r_tran ) ) then
         tmin = min(density(3),density(2))
         tmax = max(density(3),density(2))
         density_at_edge = min( tmax, max(density_at_edge,tmin) )
      end if
    end if

  end function calc_density_at_cell_edge

  !----------------------------------------------------------------------------
  !> @brief  Outputs the coefficients (a0,a1,a2) for the subgrid representation
  !!         rho(x) = a0 + a1*x + a2*x^2. Inputs are the density value for the
  !!         cell, and the left hand and right hand estimates of the density
  !!         for the cell. Given these three values a quadratic subgrid
  !!         approximation of rho can be made. Monotonicity is applied if
  !!         selected.
  !!
  !! @param[in]   density_cell_edge_left   Estimate of the density at x=0
  !! @param[in]   density_cell_edge_right  Estimate of the density at x=1
  !! @param[in]   density_of_cell          Average density of the cell
  !! @param[in]   monotone                 Monotone option to ensures no over/undershoots
  !! @param[out]  coeffs                   coeffs(1)=a0, coeffs(2)=a1, coeffs(3)=a2
  !----------------------------------------------------------------------------
  subroutine ppm_output(density_cell_edge_left,density_cell_edge_right,density_of_cell,monotone,coeffs)

    implicit none

    real(kind=r_tran),    intent(in)  :: density_cell_edge_left
    real(kind=r_tran),    intent(in)  :: density_cell_edge_right
    real(kind=r_tran),    intent(in)  :: density_of_cell
    integer(kind=i_def),  intent(in)  :: monotone
    real(kind=r_tran),    intent(out) :: coeffs(1:3)

    real(kind=r_tran) :: t1,t2,t3

    ! Calculate coefficients
    coeffs(1) = density_cell_edge_left
    coeffs(2) = -4.0_r_tran*density_cell_edge_left - 2.0_r_tran*density_cell_edge_right + 6.0_r_tran*density_of_cell
    coeffs(3) =  3.0_r_tran*density_cell_edge_left + 3.0_r_tran*density_cell_edge_right - 6.0_r_tran*density_of_cell

    ! Apply monotonicity if needed
    if ( monotone == horizontal_monotone_strict .OR. &
         monotone == vertical_monotone_strict ) then
      ! Strict monotonicity
      t1 = -0.5_r_tran*coeffs(2)/(coeffs(3) + EPS)
      if ((t1+EPS)*(1.0_r_tran+EPS-t1) > 0.0_r_tran) then
        coeffs(1) = density_of_cell
        coeffs(2) = 0.0_r_tran
        coeffs(3) = 0.0_r_tran
      end if
    else if ( monotone == horizontal_monotone_relaxed .OR. &
              monotone == vertical_monotone_relaxed ) then
      ! Relaxed monotonicity
      t1 = -0.5_r_tran*coeffs(2)/(coeffs(3) + EPS)
      if ((t1+EPS)*(1.0_r_tran+EPS-t1) > 0.0_r_tran) then
        t2 = (density_cell_edge_right-density_of_cell) * (density_of_cell-density_cell_edge_left)
        t3 = abs(density_of_cell-density_cell_edge_left) - abs(density_cell_edge_right-density_of_cell)
        if ( t2 < 0.0_r_tran ) then
          coeffs(1) = density_of_cell
          coeffs(2) = 0.0_r_tran
          coeffs(3) = 0.0_r_tran
        else
          if ( t3 < 0.0_r_tran ) then
            coeffs(1) = density_cell_edge_left
            coeffs(2) = 0.0_r_tran
            coeffs(3) = 3.0_r_tran*(density_of_cell - density_cell_edge_left)
          else
            coeffs(1) = -2.0_r_tran*density_cell_edge_right + 3.0_r_tran*density_of_cell
            coeffs(2) =  6.0_r_tran*density_cell_edge_right - 6.0_r_tran*density_of_cell
            coeffs(3) = -3.0_r_tran*density_cell_edge_right + 3.0_r_tran*density_of_cell
          end if
        end if
      end if
    end if

  end subroutine ppm_output

end module subgrid_rho_mod
