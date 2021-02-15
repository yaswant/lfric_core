! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_create_lfric_fields_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int64

! lfric modules
USE constants_mod,                 ONLY : i_def
USE field_mod,                     ONLY : field_type
USE field_parent_mod,              ONLY:  read_interface, write_interface,     &
                                          field_parent_proxy_type
use lfric_xios_read_mod,           ONLY : read_field_single_face,              &
                                          read_field_face
use lfric_xios_write_mod,          ONLY : write_field_single_face,             &
                                          write_field_face
USE finite_element_config_mod,     ONLY : element_order
USE function_space_collection_mod, ONLY : function_space_collection
USE field_collection_mod,          ONLY : field_collection_type
USE fs_continuity_mod,             ONLY : W3, Wtheta
USE log_mod,                       ONLY : LOG_LEVEL_INFO, LOG_LEVEL_ERROR,     &
                                          log_event, log_scratch_space

! lfricinp modules
USE lfricinp_stashmaster_mod,   ONLY: get_stashmaster_item, levelt,       &
                                           rho_levels, theta_levels,           &
                                           single_level, soil_levels,          &
                                           pseudt
USE lfricinp_stash_to_lfric_map_mod, ONLY: get_field_name
USE lfricinp_grid_type_mod,          ONLY: lfricinp_grid_type
USE lfricinp_um_level_codes_mod,     ONLY: lfricinp_get_num_levels,            &
                                           lfricinp_get_num_pseudo_levels

! shumlib modules
USE f_shum_file_mod, ONLY: shum_file_type

IMPLICIT NONE
PRIVATE
PUBLIC :: lfricinp_create_lfric_fields

CONTAINS

SUBROUTINE lfricinp_create_lfric_fields(mesh_id, twod_mesh_id,  &
                                        field_collection, stash_list, &
                                        um_grid, um_file)
! Description:
!  When given list of stashcodes create lfric field. Use stashcode to determine
!  which vertical level set the field is on and create appropiate function
!  space. Set read and write procedure for each field.
IMPLICIT NONE

INTEGER(KIND=i_def), INTENT(IN)                 :: mesh_id
INTEGER(KIND=i_def), INTENT(IN)                 :: twod_mesh_id
TYPE(field_collection_type), INTENT(IN OUT) :: field_collection

INTEGER(KIND=int64), INTENT(IN)  :: stash_list(:)

TYPE(lfricinp_grid_type), INTENT(IN):: um_grid
TYPE(shum_file_type), INTENT(IN) :: um_file

PROCEDURE(read_interface), POINTER  :: tmp_read_ptr => NULL()
PROCEDURE(write_interface), POINTER  :: tmp_write_ptr => NULL()

TYPE( field_type ) :: field
INTEGER(KIND=int64) :: stashcode, level_code, i
INTEGER(KIND=i_def) :: ndata

CALL log_event( 'Creating lfric fields...', LOG_LEVEL_INFO )

IF (element_order > 0) THEN
  CALL log_event( 'Finite diff fields: requires lowest order elements', &
         LOG_LEVEL_ERROR )
END IF

! Create the field collection
field_collection  =  field_collection_type(name="lfric_fields")

DO i=1, SIZE(stash_list)
  stashcode = stash_list(i)
  level_code = get_stashmaster_item(stashcode, levelt)


  ! TODO - It would be good to move this case statement into its own routine
  ! it could be used elsewhere to add fields to the collection without the
  ! do loop over stash list
  SELECT CASE (level_code)

  CASE(rho_levels) ! Stashcodes that map to W3/rho
    CALL field % initialise( vector_space =                                 &
         function_space_collection%get_fs(mesh_id, element_order, W3),      &
         name=TRIM(get_field_name(stashcode)))
    tmp_read_ptr => read_field_face
    tmp_write_ptr => write_field_face

  CASE(theta_levels) ! Stashcodes that maps to Wtheta
    CALL field % initialise( vector_space =                                 &
         function_space_collection%get_fs(mesh_id, element_order, Wtheta),  &
         name=TRIM(get_field_name(stashcode)))
    tmp_read_ptr => read_field_face
    tmp_write_ptr => write_field_face

  CASE(single_level) ! Stash that needs 2D mesh

    IF ( get_stashmaster_item(stashcode, pseudt) == 0 ) THEN
      ! Field has no pseudo levels
      ndata = 1
    ELSE
      ! Get number of pseudo levels/ndata
      ndata = lfricinp_get_num_pseudo_levels(um_grid, stashcode)
    END IF

    CALL field % initialise( vector_space =                                 &
         function_space_collection%get_fs(twod_mesh_id, element_order, W3,  &
         ndata), name=TRIM(get_field_name(stashcode)))
    ! 2D fields require a different interface
    tmp_write_ptr => write_field_single_face
    tmp_read_ptr => read_field_single_face

  CASE(soil_levels) ! Soil fields
    ! Get number of soil levels and set using the multidata argument
    ndata = lfricinp_get_num_levels(um_file, stashcode)
    CALL field % initialise( vector_space =                                 &
         function_space_collection%get_fs(twod_mesh_id, element_order, W3,  &
         ndata), name=TRIM(get_field_name(stashcode)), ndata_first=.true.)
    ! 2D fields require a different interface
    tmp_write_ptr => write_field_single_face
    tmp_read_ptr => read_field_single_face

  CASE DEFAULT
    WRITE(log_scratch_space, '(A,I0,A)')                                    &
       "Level code ", level_code, " not supported"
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
  END SELECT

  CALL field%set_read_behaviour(tmp_read_ptr)
  CALL field%set_write_behaviour(tmp_write_ptr)
  CALL log_event("Adding " // TRIM(get_field_name(stashcode)) // &
       " to field collection", LOG_LEVEL_INFO)
  CALL field_collection%add_field(field)

END DO

CALL log_event( 'Finite diff fields created', LOG_LEVEL_INFO )

END SUBROUTINE lfricinp_create_lfric_fields


END MODULE lfricinp_create_lfric_fields_mod
