!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Unifies the flux computed by two mpi regions when using an overset
!!        mesh.
!> @details The FFSL scheme computes fluxes on all edges of cells within an mpi
!!          region. These are stored in two fields (flux_x and flux_y) which
!!          correspond to the local x- and y-directions in each panel
!!          (and note these might disagree across panels, i.e. x on one panel
!!          is the same as y on a neighbouring panel).
!!          When using the extended mesh option the fluxes at the edge of each
!!          panel use an extension of that panel to compute the fluxes. This
!!          means that two neighbouring panels will each use an extension of
!!          their own panel compute the edge fluxes, i.e. the fluxes are no
!!          longer unique! This kernel fixes that issue by computing the fluxes
!!          and the edges of a panel as the average of the two estimates.

module ffsl_unify_flux_kernel_mod

  use argument_mod,       only : arg_type,              &
                                 GH_FIELD, GH_REAL,     &
                                 CELL_COLUMN, GH_WRITE, &
                                 STENCIL, CROSS2D,      &
                                 GH_READ,               &
                                 ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,      only : r_tran, i_def, r_def
  use fs_continuity_mod,  only : W2h, W2broken
  use kernel_mod,         only : kernel_type
  use log_mod,            only : LOG_LEVEL_ERROR, &
                                 log_event

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: ffsl_unify_flux_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                       &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE,  W2h),                        &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE,  W2h),                        &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,   W2broken, STENCIL(CROSS2D)), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,   W2broken, STENCIL(CROSS2D)), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,   ANY_DISCONTINUOUS_SPACE_1,   &
                                                           STENCIL(CROSS2D))  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: ffsl_unify_flux_kernel_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: ffsl_unify_flux_kernel_code

contains

  !> @brief Compute a unique FFSL flux across two differenct mpi regions
  !> @param[in]     nlayers       Number of vertical layers in the mesh
  !> @param[in,out] flux_x        Flux in the (local) x-direction
  !> @param[in,out] flux_y        Flux in the (local) y-direction
  !> @param[in]     flux_x_broken Flux in the (local) x-direction in the broken W2 space
  !> @param[in]     smap_x_size   Size of the stencil map for the flux_x_broken field
  !> @param[in]     smap_x_max    Max depth of the stencil map for the flux_x_broken field
  !> @param[in]     smap_x        Stencil dofmap for the flux_x_broken field
  !> @param[in]     flux_y_broken Flux in the (local) y-direction in the broken W2 space
  !> @param[in]     smap_y_size   Size of the stencil map for the flux_y_broken field
  !> @param[in]     smap_y_max    Max depth of the stencil map for the flux_y_broken field
  !> @param[in]     smap_y        Stencil dofmap for the flux_y_broken field
  !> @param[in]     panel_id      ID of the cubed sphere panels (1-6)
  !> @param[in]     smap_p_size   Size of the stencil map for the panel_id field
  !> @param[in]     smap_p_max    Max depth of the stencil map for the panel_id field
  !> @param[in]     smap_p        Stencil dofmap for the panel_id field
  !> @param[in]     ndf_w2h       Number of degrees of freedom for a Cell in the W2h space
  !> @param[in]     undf_w2h      Total number of degrees of freedom for the W2h space
  !> @param[in]     map_w2h       Cell dofmap for the W2h space
  !> @param[in]     ndf_w2b       Number of degrees of freedom for a Cell in the broken W2 space
  !> @param[in]     undf_w2b      Total number of degrees of freedom for the broken W2 space
  !> @param[in]     map_w2b       Cell dofmap for the broken W2 space
  !> @param[in]     ndf_id        Number of degrees of freedom for a Cell in the panel_id space
  !> @param[in]     undf_id       Total number of degrees of freedom for the panel_id space
  !> @param[in]     map_id        Cell dofmap for the panel_id space
  subroutine ffsl_unify_flux_kernel_code( nlayers,       &
                                          flux_x,        &
                                          flux_y,        &
                                          flux_x_broken, &
                                          smap_x_size,   &
                                          smap_x_max,    &
                                          smap_x,        &
                                          flux_y_broken, &
                                          smap_y_size,   &
                                          smap_y_max,    &
                                          smap_y,        &
                                          panel_id,      &
                                          smap_p_size,   &
                                          smap_p_max,    &
                                          smap_p,        &
                                          ndf_w2h,       &
                                          undf_w2h,      &
                                          map_w2h,       &
                                          ndf_w2b,       &
                                          undf_w2b,      &
                                          map_w2b,       &
                                          ndf_id,        &
                                          undf_id,       &
                                          map_id )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: undf_w2b
    integer(kind=i_def), intent(in) :: undf_id
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2b
    integer(kind=i_def), intent(in) :: ndf_id

    ! Arguments: Maps
    integer(kind=i_def), dimension(ndf_w2h), intent(in) :: map_w2h
    integer(kind=i_def), dimension(ndf_w2b), intent(in) :: map_w2b
    integer(kind=i_def), dimension(ndf_id),  intent(in) :: map_id

    ! Arguments: Stencil maps
    integer(kind=i_def),                                    intent(in) :: smap_x_max
    integer(kind=i_def), dimension(4),                      intent(in) :: smap_x_size
    integer(kind=i_def), dimension(ndf_w2b, smap_x_max, 4), intent(in) :: smap_x
    integer(kind=i_def),                                    intent(in) :: smap_y_max
    integer(kind=i_def), dimension(4),                      intent(in) :: smap_y_size
    integer(kind=i_def), dimension(ndf_w2b, smap_y_max, 4), intent(in) :: smap_y
    integer(kind=i_def),                                    intent(in) :: smap_p_max
    integer(kind=i_def), dimension(4),                      intent(in) :: smap_p_size
    integer(kind=i_def), dimension(ndf_id,  smap_p_max, 4), intent(in) :: smap_p

    ! Arguments: Fields
    real(kind=r_tran), dimension(undf_w2h), intent(inout) :: flux_x, flux_y
    real(kind=r_tran), dimension(undf_w2b), intent(in)    :: flux_x_broken, flux_y_broken
    real(kind=r_def),  dimension(undf_id),  intent(in)    :: panel_id

    ! Local indices
    integer(kind=i_def) :: k, df, owned, halo

    ! Only perform the computation if we have a full stencil (i.e. not at a Lam
    ! boundary)
    if (minval(smap_x_size) == 1) return

    ! The folling code assumes a certain layout and connectivity of the cubed
    ! sphere:
    !
    !       |----|
    !       | 5  |
    !  |----|----|----|----|
    !  | 4  | 1  | 2  | 3  |
    !  |----|----|----|----|
    !       | 6  |
    !       |----|
    !
    ! And the local edges of each panel are
    !         |------|
    !         |  4   |
    !         |1    3|
    !         |  2   |
    !  |------|------|------|------|
    !  |  3   |  4   |  4   |   3  |
    !  |4    2|1    3|1    3|4    2|
    !  |  1   |  2   |  2   |   1  |
    !  |------|------|------|------|
    !         |  1   |
    !         |2    4|
    !         |  3   |
    !         |------|


    ! Check if the owned panel is different from a neighbour cell
    ! and if so set the flux on that edge to be the average
    ! of that computed by this cell and the neighbouring cell,
    ! taking care to use the appropriate flux_x or flux_y fields
    ! where the orientation changes across a panel boundary
    owned = int(panel_id(map_id(1)),i_def)

    select case (owned)

      case(1)
        do df = 1, 4
          halo = int(panel_id(smap_p(1, 2, df)), i_def)
          select case(halo)
            case(2)
              do k = 0, nlayers-1
                flux_x(map_w2h(3)+k) = 0.5_r_tran*(flux_x_broken(map_w2b(3)+k) &
                                                 + flux_x_broken(smap_x(1,2,3)+k))
              end do
            case(4)
              do k = 0, nlayers-1
                flux_x(map_w2h(1)+k) = 0.5_r_tran*(flux_x_broken(map_w2b(1)+k) &
                                                 - flux_y_broken(smap_y(2,2,1)+k))
              end do
            case(5)
              do k = 0, nlayers-1
                flux_y(map_w2h(4)+k) = 0.5_r_tran*(flux_y_broken(map_w2b(4)+k) &
                                                 + flux_y_broken(smap_x(2,2,4)+k))
              end do
            case(6)
              do k = 0, nlayers-1
                flux_y(map_w2h(2)+k) = 0.5_r_tran*(flux_y_broken(map_w2b(2)+k) &
                                                 - flux_x_broken(smap_x(1,2,2)+k))
              end do
          end select
        end do

      case(2)
        do df = 1, 4
          halo = int(panel_id(smap_p(1, 2, df)), i_def)
          select case(halo)
            case(1)
              do k = 0, nlayers-1
                flux_x(map_w2h(1)+k) = 0.5_r_tran*(flux_x_broken(map_w2b(1)+k) &
                                                 + flux_x_broken(smap_x(3,2,1)+k))
              end do
            case(3)
              do k = 0, nlayers-1
                flux_x(map_w2h(3)+k) = 0.5_r_tran*(flux_x_broken(map_w2b(3)+k) &
                                                 - flux_y_broken(smap_y(4,2,3)+k))
              end do
            case(5)
              do k = 0, nlayers-1
                flux_y(map_w2h(4)+k) = 0.5_r_tran*(flux_y_broken(map_w2b(4)+k) &
                                                 - flux_x_broken(smap_x(3,2,4)+k))
              end do
            case(6)
              do k = 0, nlayers-1
                flux_y(map_w2h(2)+k) = 0.5_r_tran*(flux_y_broken(map_w2b(2)+k) &
                                                 + flux_y_broken(smap_y(4,2,2)+k))
              end do
          end select
        end do

      case(3)
        do df = 1, 4
          halo = int(panel_id(smap_p(1, 2, df)), i_def)
          select case(halo)
            case(2)
              do k = 0, nlayers-1
                flux_y(map_w2h(4)+k) = 0.5_r_tran*(flux_y_broken(map_w2b(4)+k) &
                                                 - flux_x_broken(smap_x(3,2,4)+k))
              end do
            case(4)
              do k = 0, nlayers-1
                flux_y(map_w2h(2)+k) = 0.5_r_tran*(flux_y_broken(map_w2b(2)+k) &
                                                 + flux_y_broken(smap_y(4,2,2)+k))
              end do
            case(5)
              do k = 0, nlayers-1
                flux_x(map_w2h(3)+k) = 0.5_r_tran*(flux_x_broken(map_w2b(3)+k) &
                                                 - flux_y_broken(smap_y(4,2,3)+k))
              end do
            case(6)
              do k = 0, nlayers-1
                flux_x(map_w2h(1)+k) = 0.5_r_tran*(flux_x_broken(map_w2b(1)+k) &
                                                 + flux_x_broken(smap_x(3,2,1)+k))
              end do
          end select
        end do

      case(4)
        do df = 1, 4
          halo = int(panel_id(smap_p(1, 2, df)), i_def)
          select case(halo)
            case(1)
              do k = 0, nlayers-1
                flux_y(map_w2h(2)+k) = 0.5_r_tran*(flux_y_broken(map_w2b(2)+k) &
                                                 - flux_x_broken(smap_x(1,2,2)+k))
              end do
            case(3)
              do k = 0, nlayers-1
                flux_y(map_w2h(4)+k) = 0.5_r_tran*(flux_y_broken(map_w2b(4)+k) &
                                                 + flux_y_broken(smap_y(2,2,4)+k))
              end do
            case(5)
              do k = 0, nlayers-1
                flux_x(map_w2h(3)+k) = 0.5_r_tran*(flux_x_broken(map_w2b(3)+k) &
                                                 + flux_x_broken(smap_x(1,2,3)+k))
              end do
            case(6)
              do k = 0, nlayers-1
                flux_x(map_w2h(1)+k) = 0.5_r_tran*(flux_x_broken(map_w2b(1)+k) &
                                                 - flux_y_broken(smap_y(2,2,1)+k))
              end do
          end select
        end do

      case(5)
        do df = 1, 4
          halo = int(panel_id(smap_p(1, 2, df)), i_def)
          select case(halo)
            case(1)
              do k = 0, nlayers-1
                flux_y(map_w2h(2)+k) = 0.5_r_tran*(flux_y_broken(map_w2b(2)+k) &
                                                 + flux_y_broken(smap_y(4,2,2)+k))
              end do
            case(2)
              do k = 0, nlayers-1
                flux_x(map_w2h(3)+k) = 0.5_r_tran*(flux_x_broken(map_w2b(3)+k) &
                                                 - flux_y_broken(smap_y(4,2,3)+k))
              end do
            case(3)
              do k = 0, nlayers-1
                flux_y(map_w2h(4)+k) = 0.5_r_tran*(flux_y_broken(map_w2b(4)+k) &
                                                 - flux_x_broken(smap_x(3,2,4)+k))
              end do
            case(4)
              do k = 0, nlayers-1
                flux_x(map_w2h(1)+k) = 0.5_r_tran*(flux_x_broken(map_w2b(1)+k) &
                                                 + flux_x_broken(smap_x(3,2,1)+k))
              end do
          end select
        end do

      case(6)
        do df = 1, 4
          halo = int(panel_id(smap_p(1, 2, df)), i_def)
          select case(halo)
            case(1)
              do k = 0, nlayers-1
                flux_x(map_w2h(1)+k) = 0.5_r_tran*(flux_x_broken(map_w2b(1)+k) &
                                                 - flux_y_broken(smap_y(2,2,1)+k))
              end do
            case(2)
              do k = 0, nlayers-1
                flux_y(map_w2h(4)+k) = 0.5_r_tran*(flux_y_broken(map_w2b(4)+k) &
                                                 + flux_y_broken(smap_y(2,2,4)+k))
              end do
            case(3)
              do k = 0, nlayers-1
                flux_x(map_w2h(3)+k) = 0.5_r_tran*(flux_x_broken(map_w2b(3)+k) &
                                                 + flux_x_broken(smap_x(1,2,3)+k))
              end do
            case(4)
              do k = 0, nlayers-1
                flux_y(map_w2h(2)+k) = 0.5_r_tran*(flux_y_broken(map_w2b(2)+k) &
                                                 - flux_x_broken(smap_x(1,2,2)+k))
              end do
          end select
        end do
     case default
      call log_event('Invalid panel',LOG_LEVEL_ERROR)
    end select

  end subroutine ffsl_unify_flux_kernel_code

end module ffsl_unify_flux_kernel_mod
