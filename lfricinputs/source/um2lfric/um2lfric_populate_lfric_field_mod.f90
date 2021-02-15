! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE um2lfric_populate_lfric_field_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : real64, int32, int64

! LFRic modules
USE field_mod, ONLY: lfric_field_type => field_type,      &
                     lfric_proxy_type => field_proxy_type

IMPLICIT NONE

PRIVATE

PUBLIC :: um2lfric_populate_lfric_field

CONTAINS

!--------------------------------------------------------------------------

SUBROUTINE um2lfric_populate_lfric_field(regridded_field, lfric_field)

USE log_mod,        ONLY: log_event, LOG_LEVEL_ERROR, LOG_LEVEL_INFO, &
                          log_scratch_space
USE extrusion_config_mod, ONLY: number_of_layers

! Description:
!  Performs the copy of data from an regridded field array of levels into
!  an LFRic field proxy data array. The LFRic field is a single 1D array
!  containing all levels
!
!  The field ordering in the LFRic array loops over each column before moving
!  onto the next horizontal point.
!
!  In the current UM2LFRic code we use a serial implementation of the LFRic
!  infrastructure. This means that there are no complications from the
!  MPI partitioning or from field stencils.
!
!  Simple example of a 2 layer mesh, with 4 points on each layer and both
!  the target and LFRic field indicies labelled for each point
!
!  Note that the target array indices repeat for each layer but the LFRic field
!  index refers to the full 3D field

!                Layer 1                                  Layer 2
!
!
!     x  regridded_1      x  regridded_3       x  regridded_1       x  regridded_3
!        LFRIC_1             LFRIC_5              LFRIC_2              LFRIC_6
!
!     x  regridded_2      x  regridded_4       x  regridded_2       x  regridded_4
!        LFRIC_3             LFRIC_7              LFRIC_4              LFRIC_8

! 2D array field containing the regridded data. First dimension corresponds
! to number of points per level and the second dimension to the number levels
REAL(KIND=real64), INTENT(IN) :: regridded_field(:,:)
! LFRic field containing a 1D array which the data needs to be copied
! across to
TYPE(lfric_field_type), INTENT(INOUT) :: lfric_field

! Data arrays
INTEGER :: level, num_levels
INTEGER :: regridded_index, len_regridded_field, lfric_index
! Proxy type needed to access the data array directly. This is the recommended
! method for such an action.
TYPE(lfric_proxy_type) :: lfric_field_proxy


INTEGER :: regridded_total_field_size
INTEGER :: lfric_total_field_size

! Create a proxy to allow access to the LFRic field's data
lfric_field_proxy = lfric_field%get_proxy()

num_levels = SIZE(regridded_field, 2)

! Get the size of a single level and multiply by total number of levels to get
! total field size
regridded_total_field_size = SIZE(regridded_field, 2) * SIZE(regridded_field, 1)
lfric_total_field_size = SIZE(lfric_field_proxy % data)

IF ( regridded_total_field_size /= lfric_total_field_size) THEN
  WRITE(log_scratch_space, '(2(A,I0))') "Mismatch between size of regridded field: " , &
       regridded_total_field_size, " and size of LFRic field: ", lfric_total_field_size
  CALL log_event( log_scratch_space, LOG_LEVEL_ERROR)
END IF

! Loop over input field layers
DO level = 1, num_levels
  len_regridded_field = SIZE(regridded_field, 1)
  ! Split up the data here and insert it into the LFRic field in the
  ! correct place. Loop over all points in the 2D regridded field that
  ! represents a single level. The first lfric array index will match the
  ! current level number
  lfric_index = level
  DO regridded_index = 1, len_regridded_field
    lfric_field_proxy%data(lfric_index) = regridded_field(regridded_index, level)
    ! Need to step by the total number of levels to get to the lfric
    ! array index that corresponds to the next regridded point
    lfric_index = lfric_index + num_levels
  END DO
END DO

END SUBROUTINE um2lfric_populate_lfric_field

!--------------------------------------------------------------------------

END MODULE um2lfric_populate_lfric_field_mod
