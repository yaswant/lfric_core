! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE um2lfric_regrid_weights_mod

USE lfricinp_um_parameters_mod, ONLY: fnamelen
USE lfricinp_regrid_weights_type_mod, ONLY: lfricinp_regrid_weights_type

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int32, int64, real64

IMPLICIT NONE

PRIVATE

! Regridding Weights
TYPE(lfricinp_regrid_weights_type), PUBLIC, TARGET :: grid_p_to_mesh_face_centre
TYPE(lfricinp_regrid_weights_type), PUBLIC, TARGET :: grid_u_to_mesh_face_centre
TYPE(lfricinp_regrid_weights_type), PUBLIC, TARGET :: grid_v_to_mesh_face_centre


PUBLIC ::  um2lfric_regrid_weightsfile_ctl, get_weights

CONTAINS

!---------------------------------------------------------

SUBROUTINE um2lfric_regrid_weightsfile_ctl()
USE um2lfric_namelist_mod, ONLY: um2lfric_config
USE lfricinp_um_grid_mod, ONLY: um_grid
IMPLICIT NONE

! P points
CALL grid_p_to_mesh_face_centre % load(                    &
     um2lfric_config%weights_file_p_to_face_centre)
CALL grid_p_to_mesh_face_centre % validate_src(            &
     INT(um_grid % num_p_points_x, KIND=int32) *      &
     INT(um_grid % num_p_points_y, KIND=int32))
CALL grid_p_to_mesh_face_centre % populate_src_address_2D( &
     INT(um_grid % num_p_points_x, KIND=int32))

! U points
CALL grid_u_to_mesh_face_centre % load(                    &
     um2lfric_config%weights_file_u_to_face_centre)
CALL grid_u_to_mesh_face_centre % validate_src(            &
     INT(um_grid % num_u_points_x, KIND=int32) *      &
     INT(um_grid % num_u_points_y, KIND=int32))
CALL grid_u_to_mesh_face_centre % populate_src_address_2D( &
     INT(um_grid % num_u_points_x, KIND=int32))

! V points
CALL grid_v_to_mesh_face_centre % load(                    &
     um2lfric_config%weights_file_v_to_face_centre)
CALL grid_v_to_mesh_face_centre % validate_src(            &
     INT(um_grid % num_v_points_x, KIND=int32) *      &
     INT(um_grid % num_v_points_y, KIND=int32))
CALL grid_v_to_mesh_face_centre % populate_src_address_2D( &
     INT(um_grid % num_v_points_x, KIND=int32))

END SUBROUTINE um2lfric_regrid_weightsfile_ctl

!---------------------------------------------------------

FUNCTION get_weights(stashcode) RESULT (weights)

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int64
! UM2LFRic modules
USE lfricinp_stashmaster_mod, ONLY: stashmaster, p_points, &
   u_points, v_points, ozone_points, land_compressed,           &
   p_points_values_over_sea
! LFRic modules
USE log_mod, ONLY: log_event, LOG_LEVEL_ERROR, log_scratch_space

IMPLICIT NONE

! Arguments
INTEGER(KIND=int64), INTENT(IN) :: stashcode
! Result
TYPE(lfricinp_regrid_weights_type), POINTER :: weights

! Local variables
INTEGER(KIND=int64) :: horiz_grid_code = 0

! Check that the STASHmaster record exists.
IF (ASSOCIATED(stashmaster(stashcode) % record)) THEN
  ! Get grid type code from STASHmaster entry
  horiz_grid_code = stashmaster(stashcode) % record % grid
ELSE
  WRITE(log_scratch_space, '(A,I0)')                             &
       "Unassociated STASHmaster record for stashcode ", stashcode
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

IF (horiz_grid_code == u_points) THEN

  weights => grid_u_to_mesh_face_centre

ELSE IF (horiz_grid_code == v_points) THEN

  weights => grid_v_to_mesh_face_centre

ELSE IF (horiz_grid_code == p_points .OR.   &
         horiz_grid_code == ozone_points .OR. &
         horiz_grid_code == land_compressed .OR. &
         horiz_grid_code == p_points_values_over_sea) THEN

   weights => grid_p_to_mesh_face_centre
ELSE
  WRITE(log_scratch_space, '(2(A,I0))') "Unsupported horizontal grid type code: ", &
       horiz_grid_code, " encountered during regrid of stashcode", stashcode
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

IF (.NOT. ALLOCATED(weights%remap_matrix)) THEN
  CALL log_event("Attempted to select unallocated weights matrix", LOG_LEVEL_ERROR)
END IF

END FUNCTION get_weights


END MODULE um2lfric_regrid_weights_mod
