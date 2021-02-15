! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE scintelapi_interface_mod
!
! This module provides the necessary routines for initialising, configuring and
! finalising the API.
!

USE log_mod,                   ONLY: log_event, log_scratch_space,             &
                                     initialise_logging, finalise_logging,     &
                                     LOG_LEVEL_INFO, LOG_LEVEL_ERROR
USE scintelapi_namelist_mod,   ONLY: scintelapi_namelist_from_cl, lfric_nl,    &
                                     scintelapi_nl, required_lfric_namelists
USE constants_def_mod,         ONLY: func_space_name_len, field_name_len,      &
                                     gen_id_len, genpar_len, field_dim_len,    &
                                     field_id_list_max_size, empty_string

IMPLICIT NONE

CONTAINS


SUBROUTINE scintelapi_initialise()
!
! This routine initialises all the necessary LFRic, XIOS, YAXT and API
! infrastructure.
!

USE field_list_mod,            ONLY: init_field_list
USE generator_library_mod,     ONLY: init_generator_lib
USE dependency_graph_list_mod, ONLY: init_dependency_graph_list
USE lfricinp_lfric_driver_mod, ONLY: lfricinp_initialise_lfric

IMPLICIT NONE

! Read namelist file names from command line
CALL scintelapi_namelist_from_cl()

! Initialise LFRic infrastructure
CALL lfricinp_initialise_lfric(                                                &
                         program_name_arg="scintelapi",                        &
                         lfric_nl_fname=lfric_nl,                              &
                         required_lfric_namelists=required_lfric_namelists     &
                              )

! Initialise the field list
CALL init_field_list()

! Initialise the generator library
CALL init_generator_lib()

! Initialise the dependency graph list
CALL init_dependency_graph_list()

END SUBROUTINE scintelapi_initialise


SUBROUTINE scintelapi_finalise()
!
! This routine finalises all the used infrastructure
!

USE lfricinp_lfric_driver_mod, ONLY: lfricinp_finalise_lfric

IMPLICIT NONE

! Finalise LFRic infrastructure.
CALL lfricinp_finalise_lfric()

END SUBROUTINE scintelapi_finalise


SUBROUTINE scintelapi_add_field(field_id, func_space, field_dim, write_name)
!
! This routine is used to add a field to the internally stored global field
! list. It also performs several checks on the input for validity and
! consistency.
!

USE finite_element_config_mod,      ONLY: element_order
USE function_space_collection_mod , ONLY: function_space_collection
USE fs_continuity_mod,              ONLY: W3, Wtheta
USE lfricinp_lfric_driver_mod,      ONLY: mesh_id, twod_mesh_id
USE field_mod,                      ONLY: field_type
USE field_list_mod,                 ONLY: no_fields, field_list,               &
                                          field_io_id_list

IMPLICIT NONE

!
! Arguments
!
! Field identifier of new field to be added to list
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: field_id

! Function space identifier of function space on which this field is defined
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: func_space

! 2D or 3D field identifier
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: field_dim

! XIOS write identifier to use as defined in iodef.xml file
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: write_name

!
! Local variables
!
! Mesh id to use
INTEGER :: m_id

! Function space to use
INTEGER :: Fspace

! Iterable
INTEGER :: l

! Logical used
LOGICAL :: l_field_id_exists, l_write_name_exists, l_field_id_present,         &
           l_write_name_present, l_field_dim_present, l_func_space_present

! User feedback
WRITE(log_scratch_space,'(A)') 'Attempt to add ' // TRIM(field_id) //          &
                               ' to global field list.'
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

!
! START input checks ...
!

! Check field id is present and not already in global field list
l_field_id_present = .false.
IF (PRESENT(field_id)) THEN
  IF (TRIM(field_id) /= empty_string) l_field_id_present = .true.
END IF

IF (l_field_id_present) THEN

  l_field_id_exists = .false.
  DO l = 1, no_fields
    IF (TRIM(field_list(l)%get_name()) == TRIM(field_id)) THEN
      l_field_id_exists = .true.
      EXIT
    END IF
  END DO

  IF (l_field_id_exists) THEN
    WRITE(log_scratch_space,'(A)') TRIM(field_id) // ' already in field list'
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
  END IF

ELSE

CALL log_event('No field id provided', LOG_LEVEL_ERROR)

END IF

! Check if any fields in global field list already has the same write name
l_write_name_present = .false.
IF (PRESENT(write_name)) THEN
  IF (TRIM(write_name) /= empty_string) l_write_name_present = .true.
END IF

IF (l_write_name_present) THEN

  l_write_name_exists = .false.
  DO l = 1, no_fields
    IF (TRIM(field_io_id_list(l)) == TRIM(write_name)) THEN
      l_write_name_exists = .true.
      EXIT
    END IF
  END DO

  IF (l_write_name_exists) THEN
    WRITE(log_scratch_space,'(A)') 'The supplied write_name ' //               &
                                   TRIM(write_name) // ' for field ' //        &
                                   TRIM(field_id) // ' already exist for ' //  &
                                   'another field in the field list'
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
  END IF

END IF

! Set the mesh id based of field dimension
l_field_dim_present = .false.
IF (PRESENT(field_dim)) THEN
  IF (TRIM(field_dim) /= empty_string) l_field_dim_present = .true.
END IF

IF (l_field_dim_present) THEN

  IF (TRIM(field_dim) == '3D') THEN
    m_id = mesh_id
  ELSE IF (TRIM(field_dim) == '2D') THEN
    m_id = twod_mesh_id
  ELSE
    WRITE(log_scratch_space,'(A,A)') "Unrecognised field dimension for field ",&
                                     field_id
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
  END IF

ELSE

CALL log_event('No field dimension provided', LOG_LEVEL_ERROR)

END IF

! Set function space based on function space id provided
l_func_space_present = .false.
IF (PRESENT(func_space)) THEN
  IF (TRIM(func_space) /= empty_string) l_func_space_present = .true.
END IF

IF (l_func_space_present) THEN

  IF (TRIM(func_space) == 'Wtheta') THEN
    Fspace = Wtheta
  ELSE IF (TRIM(func_space) == 'W3') THEN
    Fspace = W3
  ELSE
    WRITE(log_scratch_space,'(A,A)') "Unsupported function space for field ",  &
                                      field_id
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
  END IF

ELSE

CALL log_event('No function space provided', LOG_LEVEL_ERROR)

END IF

!
! Input checks DONE ...
!

! Set index of new field in the list
l = no_fields + 1

! Set new field
CALL field_list(l) % initialise(vector_space =                                 &
                                function_space_collection%get_fs(              &
                                                                m_id,          &
                                                                element_order, &
                                                                Fspace         &
                                                                ),             &
                                name = field_id)

! Set write id of the new field, if defined
IF (l_write_name_present) THEN
  field_io_id_list(l) = write_name
END IF

! Update the number of defined fields
no_fields = no_fields + 1

CALL log_event('Field successfully added', LOG_LEVEL_INFO)

END SUBROUTINE scintelapi_add_field


SUBROUTINE scintelapi_add_dependency_graph(input_fields, output_fields,        &
                                          generator, genpar)
!
! This routine is used to add a dependency graph to the internally stored global
! dependency graph list. It also performs several checks on the input for validity
! and consistency.
!

USE dependency_graph_list_mod, ONLY: no_dependency_graphs, dependency_graph_list
USE field_list_mod,            ONLY: get_field_pointer, field_list, no_fields
USE generator_library_mod,     ONLY: generator_index, generator_list,          &
                                     no_generators

IMPLICIT NONE

!
! Arguments
!
! Identifier of generator to run
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: generator

! Input field id array
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: input_fields(:)

! Output field id array
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: output_fields(:)

! Optional input parameter list to the generator
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: genpar

!
! Local variables
!
! Sizes of input and output field arrays for the dep graph
INTEGER :: no_input_fields, no_output_fields

! Iterable
INTEGER :: i, j, l

! Logical to check if input and output fields are defined in global field list
LOGICAL :: l_field_does_not_exist

! Logical to check if generator is defined in generator library
LOGICAL :: l_generator_does_not_exist

! Other logicals used in checking input
LOGICAL :: l_input_fields_present, l_output_fields_present, l_generator_present

CALL log_event('Attempting to ADD DEPENDENCY GRAPH', LOG_LEVEL_INFO)

!
! START input checks ...
!

! Check the input and output field lists are minimal, i.e. no repeated field
! ids, the fields are already in the global field list, and the output field
! list has at least one field id.
l_input_fields_present = .false.
IF (PRESENT(input_fields)) THEN
  IF (SIZE(input_fields) /= 0) l_input_fields_present = .true.
END IF

IF (l_input_fields_present) THEN
  WRITE(log_scratch_space,'(100(A,1X))')                                       &
        'INPUT FIELDS:', (TRIM(input_fields(i)), i=1,SIZE(input_fields))
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
  DO i = 1, SIZE(input_fields)

    DO j = 1, SIZE(input_fields)
      IF ( (TRIM(input_fields(j)) == TRIM(input_fields(i))) .AND. j /= i ) THEN
        WRITE(log_scratch_space, '(A)') 'Repeated field ids in input list'
        CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
      END IF
    END DO

    l_field_does_not_exist = .true.
    DO j = 1, no_fields
      IF (TRIM(field_list(j)%get_name()) == TRIM(input_fields(i))) THEN
        l_field_does_not_exist = .false.
        EXIT
      END IF
    END DO

    IF (l_field_does_not_exist) THEN
      WRITE(log_scratch_space, '(A)') 'Input field not defined'
      CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
    END IF

  END DO
END IF
!
l_output_fields_present = .false.
IF (PRESENT(output_fields)) THEN
  IF (SIZE(output_fields) /= 0) l_output_fields_present = .true.
END IF

IF (l_output_fields_present) THEN
  WRITE(log_scratch_space,'(100(A,1X))')                                       &
        'OUTPUT FIELDS:', (TRIM(output_fields(i)), i=1,SIZE(output_fields))
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
  DO i = 1, SIZE(output_fields)

    DO j = 1, SIZE(output_fields)
      IF ( (TRIM(output_fields(j)) == TRIM(output_fields(i))) .AND. j /= i ) THEN
        WRITE(log_scratch_space, '(A)') 'Repeated field ids in output list'
        CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
      END IF
    END DO

    l_field_does_not_exist = .true.
    DO j = 1, no_fields
      IF (TRIM(field_list(j)%get_name()) == TRIM(output_fields(i))) THEN
        l_field_does_not_exist = .false.
        EXIT
      END IF
    END DO

    IF (l_field_does_not_exist) THEN
      WRITE(log_scratch_space, '(A)') 'Output field not defined'
      CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
    END IF

  END DO
ELSE
  WRITE(log_scratch_space, '(A)') 'No output fields in dependency graph!'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

! Check generator id has been provided and that it exists in the generator
! library
l_generator_present = .false.
IF (PRESENT(generator)) THEN
  IF (TRIM(generator) /= empty_string) l_generator_present = .true.
END IF

IF (l_generator_present) THEN
  WRITE(log_scratch_space,'(2(A,1X))') 'GENERATOR:', TRIM(generator)
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

  l_generator_does_not_exist = .true.
  DO j = 1, no_generators
    IF (TRIM(generator_list(j)%identifier) == TRIM(generator)) THEN
      l_generator_does_not_exist = .false.
      EXIT
    END IF
  END DO

  IF (l_generator_does_not_exist) THEN
    WRITE(log_scratch_space, '(A)') 'Generator not defined!'
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
  END IF

ELSE
  WRITE(log_scratch_space, '(A)') 'No generator id provided!'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

!
! Input checks DONE ...
!

! Set index of new dep graph in the list
l = no_dependency_graphs + 1

! Set up the pointers to the input fields
IF (l_input_fields_present) THEN
  no_input_fields = SIZE(input_fields)
  ALLOCATE(dependency_graph_list(l)%input_field(no_input_fields))
  DO i = 1, no_input_fields
    dependency_graph_list(l)%input_field(i)%field_ptr =>                       &
                                            get_field_pointer(input_fields(i))
  END DO
END IF

! Set up the pointers to the output fields
no_output_fields = SIZE(output_fields)
ALLOCATE(dependency_graph_list(l)%output_field(no_output_fields))
DO i = 1, no_output_fields
  dependency_graph_list(l)%output_field(i)%field_ptr =>                        &
                                           get_field_pointer(output_fields(i))
END DO

! Set up generator for this dep graph
dependency_graph_list(l)%gen = generator_list(generator_index(generator))

! Set up generator parameter list
IF (PRESENT(genpar)) THEN
  dependency_graph_list(l)%genpar = genpar
ELSE
  dependency_graph_list(l)%genpar = empty_string
END IF

! Update the number of defined dep graphs
no_dependency_graphs = no_dependency_graphs + 1

CALL log_event('DEPENDENCY GRAPH successfully added', LOG_LEVEL_INFO)

END SUBROUTINE scintelapi_add_dependency_graph


SUBROUTINE scintelapi_add_fields_from_nl()
!
! This routine is used to add a list of fields to the internally stored global
! field list. Input is read from the scintelapi_nl namelist file as set in the
! module scintelapi_namelist_mod.
!

USE lfricinp_unit_handler_mod,          ONLY: get_free_unit

IMPLICIT NONE

! Unit number to use in namelist file reading.
INTEGER :: nml_unit

! Input string declarations used in namelists
CHARACTER(LEN=field_name_len)      :: field_id
CHARACTER(LEN=func_space_name_len) :: func_space
CHARACTER(LEN=field_dim_len)       :: field_dim
CHARACTER(LEN=field_name_len)      :: write_name
CHARACTER(LEN=1)                   :: write_to_dump

! Field definition namelist
NAMELIST /field_definitions/ field_id, func_space, field_dim, write_name,      &
                             write_to_dump

! Read namelist file for field definitions, and add said fields to internal
! field list
CALL get_free_unit(nml_unit)
OPEN(unit=nml_unit, file=TRIM(scintelapi_nl))

DO ! Loop over all field_definitions namelists

  ! Initialise field definition items
  field_id   = empty_string
  func_space = empty_string
  field_dim  = empty_string
  write_name = empty_string

  ! Read namelist items. Exit loop if EOF is reached
  READ(unit=nml_unit, nml=field_definitions, END=101)

  ! Report warning to user if write_name item and write_to_dump settings are in
  ! in conflict with each other. In case of such a conflict, default position is
  ! not to write field to dump, i.e. write_name is an empty string.
  IF ((TRIM(write_name) /= empty_string).AND.(write_to_dump == 'n')) THEN
    WRITE(log_scratch_space,'(A)') 'WARNING: Dump write name set, but '    //  &
                                   'writing to dump option turned off. '   //  &
                                   'Field will not be written to dump'
    CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
    write_name = empty_string
  END IF

  IF ((TRIM(write_name) == empty_string).AND.(write_to_dump == 'y')) THEN
    WRITE(log_scratch_space,'(A)') 'WARNING: Dump write name not set, but ' // &
                                   'writing to dump option turned on. '     // &
                                   'Field will not be written to dump'
    CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
  END IF

  ! Add field to global field list
  CALL scintelapi_add_field(field_id      = field_id,                          &
                            func_space    = func_space,                        &
                            field_dim     = field_dim,                         &
                            write_name    = write_name)

END DO ! End of loop over namelists

101 CLOSE(unit=nml_unit)

END SUBROUTINE scintelapi_add_fields_from_nl


SUBROUTINE scintelapi_add_dependency_graphs_from_nl()
!
! This routine is used to add a list of dependency graphs to the internally
! stored global dependency graph list. Input is read from the scintelapi_nl
! namelist file as set in the module scintelapi_namelist_mod.
!

USE lfricinp_unit_handler_mod,          ONLY: get_free_unit

IMPLICIT NONE

! Unit number to use in namelist file reading.
INTEGER :: nml_unit

! Iterable(s)
INTEGER :: i, k

! Allocatable arrays needed
CHARACTER(LEN=field_name_len), ALLOCATABLE :: input_fields_trimmed(:)
CHARACTER(LEN=field_name_len), ALLOCATABLE :: output_fields_trimmed(:)

! Input string declarations used in namelists
CHARACTER(LEN=field_name_len)      :: output_fields(field_id_list_max_size)
CHARACTER(LEN=field_name_len)      :: input_fields(field_id_list_max_size)
CHARACTER(LEN=gen_id_len)          :: generator
CHARACTER(LEN=genpar_len)          :: genpar

! Dependency graph definition namelist
NAMELIST /dependency_graphs/ input_fields, output_fields, generator, genpar

! Read namelist file for dependency graph definitions, and add said dependency
! graphs to internal list
CALL get_free_unit(nml_unit)
OPEN(unit=nml_unit, file=TRIM(scintelapi_nl))

DO ! Loop over all dependency_graphs namelists

  ! Initialise dependency graph items. Assume input and output field name arrays
  ! are of maximum size, which will be trimmed down later.
  input_fields  = [(empty_string, i = 1,field_id_list_max_size)]
  output_fields = [(empty_string, i = 1,field_id_list_max_size)]
  generator     = empty_string
  genpar        = empty_string

  ! Read namelist items. Exit loop if EOF is reached
  READ(unit=nml_unit, nml=dependency_graphs, END=102)

  ! Trim down input field name array to minimum size.
  !
  ! First find minimum size and allocate necessary array with that size.
  k = 0
  DO i = 1, field_id_list_max_size
    IF (TRIM(input_fields(i)) /= empty_string) k = k + 1
  END DO
  ALLOCATE(input_fields_trimmed(k))
  !
  ! Second fill allocated minimun sized array with the necessary data.
  k = 0
  DO i = 1, field_id_list_max_size
    IF (TRIM(input_fields(i)) /= empty_string) THEN
      k = k + 1
      input_fields_trimmed(k) = input_fields(i)
    END IF
  END DO

  ! Trim down output field name array to minimum size.
  !
  ! First find minimum size and allocate necessary array with that size.
  k = 0
  DO i = 1, field_id_list_max_size
    IF (TRIM(output_fields(i)) /= empty_string) k = k + 1
  END DO
  ALLOCATE(output_fields_trimmed(k))
  !
  ! Second fill allocated minimun sized array with the necessary data.
  k = 0
  DO i = 1, field_id_list_max_size
    IF (TRIM(output_fields(i)) /= empty_string) THEN
      k = k + 1
      output_fields_trimmed(k) = output_fields(i)
    END IF
  END DO

  ! Add dependendency graph to global dependency graph list
  CALL scintelapi_add_dependency_graph(input_fields  = input_fields_trimmed,   &
                                       output_fields = output_fields_trimmed,  &
                                       generator     = generator,              &
                                       genpar        = genpar)

  DEALLOCATE(input_fields_trimmed)
  DEALLOCATE(output_fields_trimmed)

END DO ! End of loop over namelists

102 CLOSE(unit=nml_unit)

END SUBROUTINE scintelapi_add_dependency_graphs_from_nl

END MODULE scintelapi_interface_mod
