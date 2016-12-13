!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief init functionality for boundary test

!> @details Handles init of prognostic and coordinate fields

module init_boundary_test_mod

  use assign_coordinate_field_mod,    only : assign_coordinate_field
  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type
  use finite_element_config_mod,      only : element_order
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W0, W1, W2, W3, Wtheta
  use boundary_test_init_fields_alg_mod, only : boundary_test_init_fields_alg
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO
  use restart_control_mod,            only : restart_type
  use mr_indices_mod,                 only : nummr

  implicit none


  contains

  subroutine init_boundary_test(mesh_id, chi, u, xi, restart)

    integer(i_def), intent(in)               :: mesh_id
    ! coordinate field
    type( field_type ), intent(inout)        :: chi(3)
    ! prognostic fields
    type( field_type ), intent(inout)        :: u, xi
    type(restart_type), intent(in)           :: restart

    integer(i_def)                           :: coord, imr

    type( field_type )                       :: theta, rho

    call log_event( 'boundary test: initialisation...', LOG_LEVEL_INFO )


    ! Calculate coordinates
    do coord = 1,3
      chi(coord) = field_type (vector_space =                                    &
                       function_space_collection%get_fs(mesh_id,element_order,W0) )
    end do
    ! Assign coordinate field
    call log_event( "boundary test: Computing W0 coordinate fields", LOG_LEVEL_INFO )
    call assign_coordinate_field(chi, mesh_id)


    ! Create prognostic fields
    xi    = field_type( vector_space = &
                       function_space_collection%get_fs(mesh_id, element_order, W1) )
    u     = field_type( vector_space = &
                       function_space_collection%get_fs(mesh_id, element_order, W2) )
    ! Create dummy fields for initial alg
    theta = field_type( vector_space = &
                       function_space_collection%get_fs(mesh_id, element_order, W0) )
    rho   = field_type( vector_space = &
                       function_space_collection%get_fs(mesh_id, element_order, W3) )

    ! Initialise prognostic fields
    call boundary_test_init_fields_alg( u, rho, theta, xi)

    call log_event( 'boundary test initialised', LOG_LEVEL_INFO )

  end subroutine init_boundary_test

end module init_boundary_test_mod
