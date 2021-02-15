! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_read_command_line_args_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY: int32

IMPLICIT NONE

PRIVATE

PUBLIC :: lfricinp_read_command_line_args

CONTAINS

SUBROUTINE lfricinp_read_command_line_args(lfricinputs_fname, lfric_fname)

USE log_mod,         ONLY: log_event, LOG_LEVEL_INFO, LOG_LEVEL_ERROR,      &
                           log_scratch_space
USE lfricinp_um_parameters_mod, ONLY: fnamelen
IMPLICIT NONE

CHARACTER(LEN=fnamelen), INTENT(OUT) :: lfricinputs_fname
CHARACTER(LEN=fnamelen), INTENT(OUT) :: lfric_fname

! Other variables
INTEGER :: arglen
INTEGER(KIND=int32) :: icode_32

! Read UM namelist filename from command line
CALL GET_COMMAND_ARGUMENT(1, lfricinputs_fname, arglen, icode_32)
SELECT CASE(icode_32)
CASE (0)
  CONTINUE
CASE (1)
  log_scratch_space = 'No lfricinputs namelist filename provided on command line'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
CASE (-1)
  log_scratch_space = &
             'The filename and path to the namelist file is too long '//     &
             'for currently compiled string declaration. Please '     //     &
             'recompile with a larger filenamelength parameter or '   //     &
             'reconsider location of the provided namelist.'
  icode_32 = 1 ! Force fatal error
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
CASE DEFAULT
  log_scratch_space = 'Unknown error reading namelist file from command line.'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END SELECT

! Read LFRic infrastructure namelist filename from command line
CALL GET_COMMAND_ARGUMENT(2, lfric_fname, arglen, icode_32)
SELECT CASE(icode_32)
CASE (0)
  CONTINUE
CASE (1)
  log_scratch_space = &
       'No LFRic infrastructure namelist filename provided on command line'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
CASE (-1)
  log_scratch_space = &
             'The filename and path to the namelist file is too long '//     &
             'for currently compiled string declaration. Please '     //     &
             'recompile with a larger filenamelength parameter or '   //     &
             'reconsider location of the provided namelist.'
  icode_32 = 1 ! Force fatal error
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
CASE DEFAULT
  log_scratch_space = 'Unknown error reading namelist file from command line.'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END SELECT
RETURN
END SUBROUTINE lfricinp_read_command_line_args


END MODULE lfricinp_read_command_line_args_mod
