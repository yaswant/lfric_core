!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Initialises functionality for transport.
!> @details Handles initialisation of wind, density and chi fields.
!>          transport_init_fields_alg() is used to initialise density.

module init_transport_mod

  use constants_mod,                  only: i_def
  use field_mod,                      only: field_type
  use field_parent_mod,               only: write_interface
  use finite_element_config_mod,      only: element_order
  use fs_continuity_mod,              only: W2, W3, Wtheta
  use runtime_constants_mod,          only: create_runtime_constants
  use function_space_mod,             only: function_space_type
  use function_space_collection_mod,  only: function_space_collection
  use io_config_mod,                  only: write_diag, &
                                            use_xios_io
  use write_methods_mod,              only: write_field_face
  use log_mod,                        only: log_event,                &
                                            LOG_LEVEL_INFO
  use transport_init_fields_alg_mod,  only: transport_init_fields_alg

  implicit none

  contains

  !> @param[in] mesh_id                Mesh-id
  !> @param[in] twod_mesh_id           2D Mesh-id
  !> @param[in,out] chi                Coordinate field
  !> @param[in] shifted_mesh_id        Mesh-id for vertically shifted coordinates
  !> @param[in,out] shifted_chi        Coordinate field for vertically shifted coordinates
  !> @param[in,out] wind_n             Wind field at timestep n
  !> @param[in,out] density            Density field
  !> @param[in,out] theta              Theta field
  !> @param[in,out] dep_pts_x          Departure points in the x-direction
  !> @param[in,out] dep_pts_y          Departure points in the y-direction
  !> @param[in,out] dep_pts_z          Departure points in the z-direction
  !> @param[in,out] increment          Density increment
  !> @param[in,out] divergence         Divergence field
  !> @param[in,out] wind_shifted       Wind field on vertically shifted W2 field
  !> @param[in,out] density_shifted    Density field on vertically shifted W3 field
  subroutine init_transport( mesh_id, twod_mesh_id, chi, shifted_mesh_id, shifted_chi, &
                             wind_n, density, theta, dep_pts_x, dep_pts_y, dep_pts_z,  &
                             increment, divergence, wind_shifted, density_shifted )

    implicit none

    integer(i_def),   intent(in)      :: mesh_id
    integer(i_def),   intent(in)      :: twod_mesh_id
    type(field_type), intent(inout)   :: chi(:)
    integer(i_def),   intent(in)      :: shifted_mesh_id
    type(field_type), intent(inout)   :: shifted_chi(:)
    type(field_type), intent(inout)   :: wind_n
    type(field_type), intent(inout)   :: density
    type(field_type), intent(inout)   :: theta
    type(field_type), intent(inout)   :: dep_pts_x
    type(field_type), intent(inout)   :: dep_pts_y
    type(field_type), intent(inout)   :: dep_pts_z
    type(field_type), intent(inout)   :: increment
    type(field_type), intent(inout)   :: divergence
    type(field_type), intent(inout)   :: wind_shifted
    type(field_type), intent(inout)   :: density_shifted

    type(function_space_type),  pointer :: function_space => null()
    procedure(write_interface), pointer :: tmp_write_ptr => null()

    call wind_n%initialise( vector_space = &
                          function_space_collection%get_fs( mesh_id, element_order, W2 ) )
    call density%initialise( vector_space = &
                          function_space_collection%get_fs( mesh_id, element_order, W3 ) )
    call theta%initialise( vector_space = &
                          function_space_collection%get_fs( mesh_id, element_order, Wtheta ) )
    call increment%initialise( vector_space = &
                          function_space_collection%get_fs( mesh_id, element_order, W3 ) )
    call divergence%initialise( vector_space = &
                          function_space_collection%get_fs( mesh_id, element_order, W3 ) )

    call dep_pts_x%initialise( vector_space = &
                          function_space_collection%get_fs( mesh_id, element_order, W2 ) )
    call dep_pts_y%initialise( vector_space = &
                          function_space_collection%get_fs( mesh_id, element_order, W2 ) )
    call dep_pts_z%initialise( vector_space = &
                          function_space_collection%get_fs( mesh_id, element_order, W2 ) )

    ! Create wind and density fields which live on the shifted coordinate field
    ! so that the finite-volume Cosmic transport scheme can be used to transport
    ! species living at Wtheta dofs.
    ! The shifted coordinate field has nlayers+1 with a half layer at the top and
    ! bottom of the column.
    call wind_shifted%initialise( vector_space = &
                          function_space_collection%get_fs( shifted_mesh_id, element_order, W2 ) )
    call density_shifted%initialise( vector_space = &
                          function_space_collection%get_fs( shifted_mesh_id, element_order, W3 ) )

    ! Create runtime_constants object.
    call create_runtime_constants( mesh_id, twod_mesh_id, chi, shifted_mesh_id, shifted_chi )

    ! Initialise density and theta fields
    call transport_init_fields_alg( wind_n,    &
                                    density,   &
                                    theta,     &
                                    dep_pts_x, &
                                    dep_pts_y, &
                                    dep_pts_z, &
                                    divergence )

    ! Set I/O behaviours for diagnostic output

    if ( write_diag .and. use_xios_io ) then
       ! Fields that are output on the XIOS face domain
       tmp_write_ptr => write_field_face
       call wind_n%set_write_behaviour( tmp_write_ptr )
       call density%set_write_behaviour( tmp_write_ptr )
       call theta%set_write_behaviour( tmp_write_ptr )
       call increment%set_write_behaviour( tmp_write_ptr )
       call divergence%set_write_behaviour( tmp_write_ptr )
    end if

    nullify( function_space, tmp_write_ptr )

  end subroutine init_transport

end module init_transport_mod
