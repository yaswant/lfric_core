! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_ancils_mod

USE constants_mod,                  ONLY : i_def, l_def
USE log_mod,                        ONLY : log_event,         &
                                           log_scratch_space, &
                                           LOG_LEVEL_INFO
USE field_collection_mod,           ONLY: field_collection_type
USE field_mod,                      ONLY : field_type
USE field_parent_mod,               ONLY : read_interface, &
                                           write_interface
USE io_config_mod,                  ONLY : use_xios_io
use lfric_xios_read_mod,            ONLY : read_field_face, &
                                           read_field_single_face
use lfric_xios_write_mod,           ONLY : write_field_face, &
                                           write_field_single_face
USE field_collection_mod,           ONLY : field_collection_type

USE function_space_mod,             ONLY : function_space_type
USE function_space_collection_mod,  ONLY : function_space_collection
USE fs_continuity_mod,              ONLY : W3

IMPLICIT NONE

LOGICAL :: l_land_area_fraction = .FALSE.
! Container for ancil fields
TYPE(field_collection_type) :: ancil_fields

PUBLIC   :: lfricinp_create_ancil_fields,         &
            lfricinp_setup_ancil_field,           &
            l_land_area_fraction

CONTAINS

! Organises fields to be read from ancils into ancil_fields collection

SUBROUTINE lfricinp_create_ancil_fields(ancil_fields, mesh_id, twod_mesh_id)

IMPLICIT NONE

TYPE( field_collection_type ), INTENT( OUT )   :: ancil_fields
INTEGER(i_def), INTENT(IN) :: mesh_id
INTEGER(i_def), INTENT(IN) :: twod_mesh_id

! Set up ancil_fields collection
WRITE(log_scratch_space,'(A,A)') "Create ancil fields: "// &
     "Setting up ancil field collection"
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
ancil_fields = field_collection_type(name='ancil_fields')

IF (l_land_area_fraction) THEN
  ! Surface ancils
  CALL log_event("Create land area fraction ancil", LOG_LEVEL_INFO)
  CALL lfricinp_setup_ancil_field("land_area_fraction", ancil_fields, mesh_id, &
       twod_mesh_id, twod=.TRUE.)
END IF

END SUBROUTINE lfricinp_create_ancil_fields

!------------------------------------------------------------

! Creates fields to be read into from ancillary files

SUBROUTINE lfricinp_setup_ancil_field( name, ancil_fields, mesh_id, &
                                       twod_mesh_id, twod, ndata )

IMPLICIT NONE

CHARACTER(*), INTENT(IN)                    :: name
TYPE(field_collection_type), INTENT(IN OUT) :: ancil_fields
INTEGER(i_def), INTENT(IN)                  :: mesh_id
INTEGER(i_def), INTENT(IN)                  :: twod_mesh_id
LOGICAL(l_def), OPTIONAL, INTENT(IN)        :: twod
INTEGER(i_def), OPTIONAL, INTENT(IN)        :: ndata

! Local variables
TYPE(field_type)           :: new_field
INTEGER(i_def)             :: ndat
INTEGER(i_def), PARAMETER  :: fs_order = 0

! Pointers
TYPE(function_space_type),       POINTER  :: w3_space => NULL()
TYPE(function_space_type),       POINTER  :: twod_space => NULL()
PROCEDURE(read_interface),       POINTER  :: tmp_read_ptr => NULL()
PROCEDURE(write_interface),      POINTER  :: tmp_write_ptr => NULL()
TYPE(field_type),                POINTER  :: tgt_ptr => NULL()

! Set field ndata if argument is present, else leave as default value
IF (PRESENT(ndata)) THEN
  ndat = ndata
ELSE
  ndat = 1
END IF

! Set up function spaces for field initialisation
w3_space   => function_space_collection%get_fs( mesh_id, fs_order, &
              W3, ndat )
twod_space => function_space_collection%get_fs( twod_mesh_id, fs_order, &
               W3, ndat )

! Create field
WRITE(log_scratch_space,'(3A,I6)') &
     "Creating new field for ", TRIM(name)
CALL log_event(log_scratch_space,LOG_LEVEL_INFO)
IF (PRESENT(twod)) THEN
  CALL new_field%initialise( twod_space, name=TRIM(name) )
ELSE
  CALL new_field%initialise( w3_space, name=TRIM(name) )
END IF

! Add the new field to the field collection
CALL ancil_fields%add_field(new_field)

! Get a field pointer from the collection
tgt_ptr => ancil_fields%get_field(name)

! Set up field read behaviour for 2D and 3D fields
IF (PRESENT(twod)) THEN
  tmp_read_ptr => read_field_single_face
  tmp_write_ptr => write_field_single_face
ELSE
  tmp_read_ptr => read_field_face
  tmp_write_ptr => write_field_face
END IF

! Set field read behaviour for target field
CALL tgt_ptr%set_read_behaviour(tmp_read_ptr)
CALL tgt_ptr%set_write_behaviour(tmp_write_ptr)


! Nullify pointers
NULLIFY(w3_space)
NULLIFY(twod_space)
NULLIFY(tmp_read_ptr)
NULLIFY(tmp_write_ptr)
NULLIFY(tgt_ptr)

END SUBROUTINE lfricinp_setup_ancil_field

END MODULE lfricinp_ancils_mod
