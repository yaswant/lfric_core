!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Create LBC fields.
!> @details Create LBC field collection and add fields.
module create_lbcs_mod

  use constants_mod,              only : i_def, l_def
  use log_mod,                    only : log_event,                  &
                                         LOG_LEVEL_INFO
  use field_mod,                  only : field_type
  use field_collection_mod,       only : field_collection_type
  use field_parent_mod,           only : write_interface,            &
                                         checkpoint_write_interface, &
                                         checkpoint_read_interface
  use finite_element_config_mod,  only : element_order
  use fs_continuity_mod,          only : W0, W2, W3, Wtheta
  use function_space_collection_mod, &
                                  only : function_space_collection
  use function_space_mod,         only : function_space_type
  use pure_abstract_field_mod,    only : pure_abstract_field_type

  implicit none

  private :: add_lbc_field
  public  :: create_lbc_fields

  contains

  !> @brief   Create and add LBC fields.
  !> @details Create the lateral boundary condition field collection.
  !!          On every timestep these fields will be updated and used by the
  !!          limited area model.
  !> @param[in]     mesh_id           The current 3d mesh identifier.
  !> @param[in,out] depository        Main collection of all fields in memory.
  !> @param[in,out] prognostic_fields The prognostic variables in the model.
  !> @param[in,out] lbc_fields        The lbc_fields used on every timestep.
  subroutine create_lbc_fields( mesh_id, depository, prognostic_fields, &
                                lbc_fields )

    implicit none

    integer(i_def),              intent(in)    :: mesh_id
    type(field_collection_type), intent(inout) :: depository
    type(field_collection_type), intent(inout) :: prognostic_fields
    type(field_collection_type), intent(inout) :: lbc_fields
    type(function_space_type),   pointer       :: wtheta_space => null(), &
                                                  w3_space => null(),     &
                                                  w2_space => null()
    logical(l_def)                             :: checkpoint_restart_flag

    call log_event( 'Create LBC fields', LOG_LEVEL_INFO )

    lbc_fields = field_collection_type( name='lbc_fields' )

    wtheta_space => function_space_collection%get_fs( &
                                              mesh_id, element_order, Wtheta )
    w3_space     => function_space_collection%get_fs( &
                                              mesh_id, element_order, W3 )
    w2_space     => function_space_collection%get_fs( &
                                              mesh_id, element_order, W2 )

    ! Checkpoint all necessary LBC fields.
    ! At the moment, LBC fields will be set analytically and so all fields
    ! need to be checkpointed.
    ! When we add an option to read in LBC fields from file,  some
    ! LBC fields will no longer need to be checkpointed and this flag can be
    ! set appropriately for each field.
    checkpoint_restart_flag = .true.

    call add_lbc_field( lbc_fields, depository, prognostic_fields, &
         "lbc_theta", wtheta_space, checkpoint_restart_flag )

    call add_lbc_field( lbc_fields, depository, prognostic_fields, &
         "lbc_u", w2_space, checkpoint_restart_flag )

    call add_lbc_field( lbc_fields, depository, prognostic_fields, &
         "lbc_rho", w3_space, checkpoint_restart_flag )

    call add_lbc_field( lbc_fields, depository, prognostic_fields, &
         "lbc_exner", w3_space, checkpoint_restart_flag )

    call add_lbc_field( lbc_fields, depository, prognostic_fields, &
         "boundary_u_diff", w2_space, checkpoint_restart_flag )

    call add_lbc_field( lbc_fields, depository, prognostic_fields, &
         "boundary_u_driving", w2_space, checkpoint_restart_flag )

  end subroutine create_lbc_fields

  !> @brief Add an LBC field to the depository, add pointers to the LBC
  !!        field collection and prognostic fields, and set its write,
  !!        checkpoint and restart behaviour.
  !> @param[in,out] lbc_fields               LBC fields field collection
  !> @param[in,out] depository               Depository field collection
  !> @param[in,out] prognostic_fields        Prognostic field collection
  !> @param[in]     name                     Name of the LBC field to be added
  !> @param[in]     vector_space             Function space of the field
  !> @param[in]     checkpoint_restart_flag  Flag to set checkpoint behaviour
  subroutine add_lbc_field( lbc_fields, depository, prognostic_fields, &
                            name, vector_space, checkpoint_restart_flag )

    use io_config_mod,        only : use_xios_io,             &
                                     write_diag
    use lfric_xios_read_mod,  only : checkpoint_read_xios
    use lfric_xios_write_mod, only : write_field_node, &
                                     write_field_face, &
                                     checkpoint_write_xios
    use io_mod,               only : checkpoint_write_netcdf, &
                                     checkpoint_read_netcdf

    implicit none

    type(field_collection_type),        intent(inout) :: lbc_fields
    type(field_collection_type),        intent(inout) :: depository
    type(field_collection_type),        intent(inout) :: prognostic_fields
    character(*),                       intent(in)    :: name
    type(function_space_type), pointer, intent(in)    :: vector_space
    logical(l_def),                     intent(in)    :: checkpoint_restart_flag
    type(field_type)                                  :: new_field
    class(pure_abstract_field_type),    pointer       :: field_ptr => null()

    ! pointers for xios write interface
    procedure(write_interface), pointer :: write_diag_behaviour => null()
    procedure(checkpoint_write_interface), &
                                pointer :: checkpoint_write_behaviour => null()
    procedure(checkpoint_read_interface), &
                                pointer :: checkpoint_read_behaviour => null()

    call new_field%initialise( vector_space, name=trim(name) )

    ! Set diagnostic write behaviour.
    if ( use_xios_io .and. write_diag ) then

      if ( new_field%which_function_space() == W0 ) then
        write_diag_behaviour => write_field_node
      else
        write_diag_behaviour => write_field_face
      endif

      call new_field%set_write_behaviour( write_diag_behaviour )

    endif

    ! Set checkpoint write and read behaviour.
    if ( checkpoint_restart_flag ) then
      if ( use_xios_io ) then
        checkpoint_write_behaviour => checkpoint_write_xios
        checkpoint_read_behaviour  => checkpoint_read_xios
      else
        checkpoint_write_behaviour => checkpoint_write_netcdf
        checkpoint_read_behaviour  => checkpoint_read_netcdf
      endif

      call new_field%set_checkpoint_write_behaviour(checkpoint_write_behaviour)
      call new_field%set_checkpoint_read_behaviour(checkpoint_read_behaviour)
    endif

    ! Add field to depository, and add references to field collections.
    call depository%add_field( new_field )

    field_ptr => depository%get_field( name )

    call lbc_fields%add_reference_to_field( field_ptr )

    if ( checkpoint_restart_flag ) then
      call prognostic_fields%add_reference_to_field( field_ptr )
    endif

  end subroutine add_lbc_field

end module create_lbcs_mod
