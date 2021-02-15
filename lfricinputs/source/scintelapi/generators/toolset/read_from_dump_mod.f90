! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE read_from_dump_mod
!
! This module contains a generator to read a field from dump.
!

USE dependency_graph_mod, ONLY: dependency_graph

IMPLICIT NONE

CONTAINS

SUBROUTINE read_from_dump(dep_graph)
!
! This generator reads the output field in the dependency graph from dump
!

USE gen_io_check_mod,       ONLY: gen_io_check
use lfric_xios_read_mod,    ONLY: read_field_face, read_field_single_face
USE field_mod,              ONLY: field_proxy_type
USE field_parent_mod,       ONLY: read_interface
USE log_mod,                ONLY: log_event, log_scratch_space, LOG_LEVEL_ERROR
USE constants_def_mod,      ONLY: genpar_len, field_name_len

IMPLICIT NONE

!
! Argument definitions:
!
! Dependency graph to be processed
CLASS(dependency_graph), INTENT(IN OUT) :: dep_graph

!
! Local variables
!
! IO procedure pointers
PROCEDURE(read_interface), POINTER :: tmp_read_ptr
!
! Number of layers in field
INTEGER :: no_layers
!
! Field proxies to use
TYPE(field_proxy_type) :: field_proxy
!
! Parameter list
CHARACTER(LEN=genpar_len) :: parlist
!
! Id of field to read from dump
CHARACTER(LEN=field_name_len) field_dump_id

! Error code for reading parameter list
INTEGER :: ioerr

!
! Perform some initial input checks
!
CALL gen_io_check(                                                             &
                  dep_graph=dep_graph,                                         &
                  input_field_no=0,                                            &
                  output_field_no=1,                                           &
                  parameter_no=1                                               &
                 )
!
! Done with initial input checks
!

! Get id of field to read from dump
parlist = dep_graph % genpar
READ(parlist,'(A)',IOSTAT=ioerr) field_dump_id
IF (ioerr /= 0 ) THEN
  WRITE(log_scratch_space,'(A)') 'Error occured when parsing parameter ' //    &
                                 'list in routine READ_FROM_DUMP. ' //         &
                                 'Input should be a single character string.'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

! Set read pointer based on whether field is 2D or 3D
field_proxy = dep_graph % output_field(1) % field_ptr % get_proxy()
no_layers = field_proxy % vspace % get_nlayers()

! Set up field dump ID with "read_" prefix for XIOS
field_dump_id = "read_" // TRIM(field_dump_id)

IF (no_layers == 1) THEN   ! For a 2D field ...
  tmp_read_ptr => read_field_single_face

ELSE                       ! For a 3D field ...
  tmp_read_ptr => read_field_face

END IF

! Read field from dump
CALL dep_graph % output_field(1) % field_ptr % set_read_behaviour(tmp_read_ptr)
CALL dep_graph % output_field(1) % field_ptr % read_field(TRIM(field_dump_id))

! Nullify read pointer
NULLIFY(tmp_read_ptr)

END SUBROUTINE read_from_dump

END MODULE read_from_dump_mod
