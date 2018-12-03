!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the UM convection scheme.
!>
module conv_kernel_mod

  use argument_mod,           only : arg_type,                       &
                                     GH_FIELD, GH_READ, GH_WRITE,    &
                                     CELLS
  use constants_mod,          only : i_def, i_um, r_def, r_um
  use fs_continuity_mod,      only : W3, Wtheta
  use kernel_mod,             only : kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: conv_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                &
         arg_type(GH_FIELD, GH_WRITE, WTHETA),         &
         arg_type(GH_FIELD, GH_WRITE, WTHETA),         &
         arg_type(GH_FIELD, GH_READ,  WTHETA),         &
         arg_type(GH_FIELD, GH_READ,  WTHETA),         &
         arg_type(GH_FIELD, GH_READ,  WTHETA),         &
         arg_type(GH_FIELD, GH_READ,  W3)              &
        /)
    integer :: iterates_over = CELLS
  contains
    procedure, nopass ::conv_code
  end type conv_kernel_type

  !-----------------------------------------------------------------------------
  ! Constructors
  !-----------------------------------------------------------------------------

  ! Overload the default structure constructor for function space
  interface conv_kernel_type
    module procedure conv_kernel_constructor
  end interface

  public conv_code

contains

  type(conv_kernel_type) function conv_kernel_constructor() result(self)
    return
  end function conv_kernel_constructor

  !> @brief Interface to the Lambert-Lewis convection scheme
  !> @details The Lambert-Lewis convection scheme is a simple
  !>           convection parametrization that mixes theta and moisture
  !>           as documented in UMDP41
  !! @param[in]  nlayers      Number of layers
  !! @param[out] dt_conv      Convection temperature increment
  !! @param[out] dmv_conv     Convection vapour increment
  !! @param[in]  theta_star   Potential temperature predictor after advection
  !! @param[in]  m_v          Vapour mixing ration after advection
  !! @param[in]  exner_in_wth Exner pressure field in wth space
  !! @param[in]  exner_in_w3  Exner pressure field in density space
  !! @param[in]  ndf_wth      Number of degrees of freedom per cell for potential temperature space
  !! @param[in]  undf_wth     Number unique of degrees of freedom  for potential temperature space
  !! @param[in]  map_wth      Dofmap for the cell at the base of the column for potential temperature space
  !! @param[in]  ndf_w3       Number of degrees of freedom per cell for density space
  !! @param[in]  undf_w3      Number unique of degrees of freedom  for density space
  !! @param[in]  map_w3       Dofmap for the cell at the base of the column for density space
  subroutine conv_code(nlayers,      &
                       dt_conv,      &
                       dmv_conv,     &
                       theta_star,   &
                       m_v,          &
                       exner_in_wth, &
                       exner_in_w3,  &
                       ndf_wth,      &
                       undf_wth,     &
                       map_wth,      &
                       ndf_w3,       &
                       undf_w3,      &
                       map_w3)

    !---------------------------------------
    ! UM modules
    !---------------------------------------
    use llcs, only: llcs_control
    use nlsizes_namelist_mod, only: row_length, rows
    use planet_constants_mod, only: p_zero, kappa

    implicit none
    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, ndf_w3
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3

    integer(kind=i_def), dimension(ndf_wth), intent(in) :: map_wth
    integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3

    real(kind=r_def), dimension(undf_wth), intent(out)  :: dt_conv, dmv_conv

    real(kind=r_def), dimension(undf_wth), intent(in)   :: theta_star, &
                                                           m_v,        &
                                                           exner_in_wth

    real(kind=r_def), dimension(undf_w3),  intent(in)   :: exner_in_w3

    ! Local variables for the kernel
    integer(kind=i_def) :: k

    real(r_um), dimension(row_length,rows,nlayers) :: theta_conv, q_conv, &
         theta_inc, q_inc, qcl_inc, cf_liquid_inc, p_layer_centres
    real(r_um), dimension(row_length,rows,nlayers+1) :: p_layer_boundaries
    real(r_um), dimension(row_length,rows) ::  conv_rain, p_star

    !-----------------------------------------------------------------------
    ! Initialise variables required from input fields
    !-----------------------------------------------------------------------
    do k = 1, nlayers
      theta_conv(1,1,k) = theta_star(map_wth(1) + k)
      q_conv(1,1,k) = m_v(map_wth(1) + k)
      p_layer_boundaries(1,1,k) = p_zero*(exner_in_w3(map_w3(1) + k-1))**(1.0_r_def/kappa)
      p_layer_centres(1,1,k) = p_zero*(exner_in_wth(map_wth(1) + k))**(1.0_r_def/kappa)
    end do
    ! Over-write p_layer_boundaries at lowest level with surface pressure
    p_layer_boundaries(1,1,1) = p_zero*(exner_in_wth(map_wth(1)))**(1.0_r_def/kappa)
    p_star(1,1) = p_layer_boundaries(1,1,1)
    ! Initialised p_layer_boundaries at top of atmosphere
    p_layer_boundaries(1,1,nlayers+1) = 0.0_r_um
    ! Initialise increments and output fields to zero
    theta_inc(:,:,:) = 0.0_r_um
    q_inc(:,:,:) = 0.0_r_um
    qcl_inc(:,:,:) = 0.0_r_um
    cf_liquid_inc(:,:,:) = 0.0_r_um
    conv_rain(:,:) = 0.0_r_um

    !-----------------------------------------------------------------------
    ! Call the convection scheme
    !-----------------------------------------------------------------------
    call llcs_control(theta_conv, q_conv, p_star, p_layer_boundaries, &
         p_layer_centres, theta_inc, q_inc, qcl_inc, cf_liquid_inc, conv_rain)

    !-----------------------------------------------------------------------
    ! Update fields to pass out
    !-----------------------------------------------------------------------
    do k = 1, nlayers
      ! "Increments" passed out are actually updated full fields, so
      ! calculate the real increment
      theta_inc(1,1,k) = theta_inc(1,1,k) - theta_conv(1,1,k)
      q_inc(1,1,k) = q_inc(1,1,k) - q_conv(1,1,k)
      ! Increments to pass out
      dt_conv(map_wth(1) + k) = theta_inc(1,1,k)*exner_in_wth(map_wth(1) + k)
      dmv_conv(map_wth(1) + k) = q_inc(1,1,k)
    end do
    ! Set lowest level value
    dt_conv(map_wth(1)) = dt_conv(map_wth(1) + 1)
    dmv_conv(map_wth(1)) = dmv_conv(map_wth(1) + 1)

  end subroutine conv_code

end module conv_kernel_mod
