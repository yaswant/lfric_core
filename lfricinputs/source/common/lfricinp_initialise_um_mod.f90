! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! This module has an underlying dependency on LFRic's logging module

MODULE lfricinp_initialise_um_mod

! Intrinsic modules
USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_BOOL

! Shumlib modules
USE f_shum_file_mod, ONLY: shum_file_type
USE f_shum_field_mod, ONLY: shum_field_type

! UM2LFRic modules
USE lfricinp_check_shumlib_status_mod, ONLY: shumlib
USE lfricinp_um_parameters_mod, ONLY: fnamelen, um_integer64

! LFRic modules
USE log_mod,     ONLY: log_event, LOG_LEVEL_INFO, LOG_LEVEL_ERROR, &
                       log_scratch_space

IMPLICIT NONE

PRIVATE

PUBLIC :: lfricinp_initialise_um, find_field_by_stashcode, &
     lfricinp_finalise_um, um_input_file

TYPE(shum_file_type), SAVE :: um_input_file

CONTAINS

! DEPENDS ON: c_shum_byteswap.o
! This is required to force fcm-make to compile the C code; whilst the built-in
! dependency analyser successfully works out that it needs to compile the
! Fortran side of the byte-swapping code, it requires an explicit statement
! to force it to compile the C part of the byte-swapping code. This is
! currently the approved way of linking Fortran and C in fcm-make.

!-------------------------------------------------------------------------------

SUBROUTINE lfricinp_initialise_um(fname)

  IMPLICIT NONE

  CHARACTER(LEN=fnamelen), INTENT(IN) :: fname
  CHARACTER(LEN=*), PARAMETER :: routinename='lfricinp_initialise_um'

  ! Load the UM file
  CALL log_event('Loading file '//TRIM(fname), LOG_LEVEL_INFO)
  CALL shumlib(routinename//'::open_file', um_input_file%open_file(fname),&
                print_on_success=.TRUE._C_BOOL)

  CALL shumlib(routinename//'::open_file', um_input_file%read_header(),  &
                print_on_success=.TRUE._C_BOOL)

END SUBROUTINE lfricinp_initialise_um

!-------------------------------------------------------------------------------

FUNCTION find_field_by_stashcode(stashcode_in)
  ! This returns the _first_ field in the input file which matches a given
  ! stash code, and aborts if if no matching fields are found.
  USE f_shum_lookup_indices_mod, ONLY: lbuser4
  USE lfricinp_um_parameters_mod, ONLY: um_imdi, msglen

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: stashcode_in
  INTEGER(KIND=um_integer64) :: stashcode
  CHARACTER(LEN=16)          :: timestring
  TYPE(shum_field_type), ALLOCATABLE :: found_fields(:)
  TYPE(shum_field_type) :: find_field_by_stashcode
  CHARACTER(LEN=*), PARAMETER :: routinename='find_field_by_stashcode'
  stashcode = stashcode_in
  IF (ALLOCATED(found_fields)) THEN
    DEALLOCATE(found_fields)
  END IF

  CALL shumlib(routinename//'::find_fields_in_file',   &
       um_input_file%find_fields_in_file(found_fields, &
       max_returned_fields=1_um_integer64,             &
       stashcode=stashcode),                           &
       print_on_success=.TRUE._C_BOOL)

  find_field_by_stashcode = found_fields(1)
  CALL shumlib(routinename//'::find_fields_in_file', &
       find_field_by_stashcode%get_timestring(timestring))

  WRITE(log_scratch_space, '(A,A)') 'Chosen field has validity time: ', &
       timestring
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

END FUNCTION find_field_by_stashcode

!-------------------------------------------------------------------------------
SUBROUTINE lfricinp_finalise_um()
  IMPLICIT NONE
  CHARACTER(LEN=*), PARAMETER :: routinename='lfricinp_finalise_um'

  ! Close file
  CALL shumlib(routinename//'::close_file', um_input_file%close_file() )

END SUBROUTINE lfricinp_finalise_um
!-------------------------------------------------------------------------------
END MODULE lfricinp_initialise_um_mod
