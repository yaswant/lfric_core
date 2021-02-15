! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE um2lfric_main_loop_mod


IMPLICIT NONE

PRIVATE

PUBLIC :: um2lfric_main_loop

CONTAINS

SUBROUTINE um2lfric_main_loop()

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : real64, int32, int64

! lfricinputs modules
USE lfricinp_check_shumlib_status_mod, ONLY: shumlib
USE lfricinp_regrid_weights_type_mod, ONLY: lfricinp_regrid_weights_type
USE lfricinp_stash_to_lfric_map_mod, ONLY: get_field_name
USE lfricinp_lfric_driver_mod, ONLY: lfric_fields

! um2lfric modules
USE lfricinp_initialise_um_mod, ONLY: um_input_file
USE um2lfric_regrid_weights_mod, ONLY: get_weights
USE um2lfric_namelist_mod, ONLY: um2lfric_config
USE um2lfric_post_process_fields_mod, ONLY: um2lfric_post_process_fields
USE um2lfric_populate_lfric_field_mod, ONLY: um2lfric_populate_lfric_field
USE um2lfric_apply_masked_field_adjustments_mod, ONLY: &
     um2lfric_apply_masked_field_adjustments

! shumlib modules
USE f_shum_field_mod, ONLY: shum_field_type

! lfric modules
USE field_mod, ONLY: lfric_field_type => field_type
USE mesh_mod,        ONLY: mesh_type
USE partition_mod,   ONLY: partition_type
USE time_config_mod, ONLY: timestep_end
USE log_mod,         ONLY: log_event, LOG_LEVEL_INFO, log_scratch_space
IMPLICIT NONE

! Array of shumlib field objects that will be returned from UM file
TYPE(shum_field_type), ALLOCATABLE  :: um_input_fields(:)

! Intermediate target field for regridding
REAL(KIND=real64), ALLOCATABLE :: regridded_field(:,:)

! Pointers to lfric objects
TYPE(lfric_field_type), POINTER :: lfric_field
TYPE(mesh_type), POINTER :: lfric_mesh
TYPE(lfricinp_regrid_weights_type), POINTER :: weights

! Other variables
INTEGER(KIND=int64) :: stashcode
! Iterators
INTEGER :: i_field, level
! Number levels of UM field
INTEGER(KIND=int32) :: num_levels

!-------------------------------------------------------------------------------
! Main loop over fields
!-------------------------------------------------------------------------------
WRITE(log_scratch_space, '(A,I0,A)') 'Will process ', &
     um2lfric_config%num_fields, ' fields'
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
DO i_field = 1, um2lfric_config%num_fields
  stashcode = um2lfric_config%stash_list(i_field)
  WRITE(log_scratch_space, '(A,I0)') 'Processing STASH code: ', stashcode
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
  WRITE(log_scratch_space,'(A,A)') 'LFRic field name: ',TRIM(get_field_name(stashcode))
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

  !-----------------------------------------------------------------------------
  ! Setup source field, get weights, get lfric field and allocate regridded
  ! field array
  !-----------------------------------------------------------------------------
  ! Get UM field array
  CALL shumlib("um2lfric::find_fields_in_file",                    &
                um_input_file%find_fields_in_file(um_input_fields, &
                stashcode = stashcode,      &
                lbproc = 0_int64))

  weights => get_weights(stashcode)
  lfric_field => lfric_fields%get_field(get_field_name(stashcode))
  lfric_mesh => lfric_field % get_mesh()

  num_levels = SIZE(um_input_fields)
  ALLOCATE(regridded_field(lfric_mesh%get_ncells_2d(), num_levels))
  DO level = 1, num_levels
    !-----------------------------------------------------------------------------
    ! Perform regridding
    !-----------------------------------------------------------------------------
    CALL weights % validate_src(SIZE(um_input_fields(level)%rdata))
    CALL weights % regrid_src_2d_dst_1d (src=um_input_fields(level)%rdata, &
                           dst=regridded_field(:, level))
    !-----------------------------------------------------------------------------
    ! Perform post masked adjustments
    !-----------------------------------------------------------------------------
    CALL um2lfric_apply_masked_field_adjustments(stashcode, &
                                                 src=um_input_fields(level)%rdata, &
                                                 dst=regridded_field(:, level))
  END DO ! loop over levels

  ! Tidy up input field memory
  IF (ALLOCATED(um_input_fields)) DEALLOCATE(um_input_fields)

  !-----------------------------------------------------------------------------
  ! Copy the regridded data to the lfric field
  !-----------------------------------------------------------------------------
  CALL um2lfric_post_process_fields(regridded_field, stashcode)

  CALL um2lfric_populate_lfric_field(regridded_field, lfric_field)
  DEALLOCATE(regridded_field)

  !-----------------------------------------------------------------------------
  ! Write fields to XIOS restart dump
  !-----------------------------------------------------------------------------
  WRITE(log_scratch_space,'(A)') "Checkpointing " // &
       TRIM(get_field_name(stashcode))
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
  CALL lfric_field%write_field(get_field_name(stashcode))

  weights => NULL()
  lfric_field => NULL()
  lfric_mesh => NULL()

END DO ! loop over stashcodes

END SUBROUTINE um2lfric_main_loop

END MODULE um2lfric_main_loop_mod
