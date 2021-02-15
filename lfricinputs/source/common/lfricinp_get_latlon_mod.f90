! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_get_latlon_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int32, int64, real64

USE log_mod, ONLY: log_scratch_space, log_event, LOG_LEVEL_ERROR

IMPLICIT NONE

PRIVATE

PUBLIC :: get_um_grid_coords, get_lfric_mesh_coords

CONTAINS

SUBROUTINE get_um_grid_coords(grid_type, idx, idy, lon, lat)
!
! This routine returns the latitude and longitude of a UM grid point with
! indices IDX and IDY. Note this routine is not compatible with variable
! resolution LAMs!
!
USE constants_mod,        ONLY: degrees_to_radians
USE lfricinp_um_grid_mod, ONLY: um_grid

IMPLICIT NONE

!
! Argument(s)
!
CHARACTER(LEN=*),    INTENT(IN)  :: grid_type
INTEGER(KIND=int64), INTENT(IN)  :: idx, idy
REAL(KIND=real64),   INTENT(OUT) :: lat, lon

IF (grid_type == 'p' ) THEN

  lon = um_grid%p_origin_x + REAL(idx-1) * um_grid%spacing_x
  lat = um_grid%p_origin_y + REAL(idy-1) * um_grid%spacing_y

ELSE IF (grid_type == 'u') THEN

  lon = um_grid%u_origin_x + REAL(idx-1) * um_grid%spacing_x
  lat = um_grid%u_origin_y + REAL(idy-1) * um_grid%spacing_y

ELSE IF (grid_type == 'v') THEN

  lon = um_grid%v_origin_x + REAL(idx-1) * um_grid%spacing_x
  lat = um_grid%v_origin_y + REAL(idy-1) * um_grid%spacing_y

ELSE

  WRITE(log_scratch_space,'(A)') 'UM grid code ' // grid_type // ' was not ' //&
                                 'recognised in routine GET_UM_GRID_COORDS'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)

END IF

! Readjust longitude range from [0,360] to [-180,180]
IF (lon > 180.0) lon = lon - 360.0

! Convert longitude and latitude to radians
lon = lon * degrees_to_radians
lat = lat * degrees_to_radians

END SUBROUTINE get_um_grid_coords


SUBROUTINE get_lfric_mesh_coords(cell_lid, lon, lat)
!
! This routine returns the latitude and longitude at the centre of a LFRic mesh
! cell with local cell id CELL_LID
!
USE mesh_collection_mod,        ONLY: mesh_collection
USE mesh_mod,                   ONLY: mesh_type
USE global_mesh_collection_mod, ONLY: global_mesh_collection
USE global_mesh_mod,            ONLY: global_mesh_type
USE lfricinp_lfric_driver_mod,  ONLY: mesh_id

IMPLICIT NONE

!
! Argument(s)
!
INTEGER,           INTENT(IN)  :: cell_lid
REAL(KIND=real64), INTENT(OUT) :: lon, lat

! Local variables
TYPE(mesh_type),        POINTER  :: mesh => NULL()
TYPE(global_mesh_type), POINTER  :: global_mesh => NULL()
INTEGER(KIND=int32), ALLOCATABLE :: verts(:)
INTEGER(KIND=int32)              :: i
REAL(KIND=real64)                :: vert_coords(2)

mesh => mesh_collection%get_mesh(mesh_id)
global_mesh => global_mesh_collection%get_global_mesh(mesh%get_global_mesh_id())
ALLOCATE(verts(global_mesh%get_nverts_per_cell()))
CALL global_mesh%get_vert_on_cell(mesh%get_cell_gid(cell_lid), verts)

!
! NOTE: Below we take the cell centre coordinates to be the average of the
! vertices' coordinates making up the cell. This is due to the fact the LFRic
! infrastucture does not provide a means of accessing the cell centre
! coordinates from the global mesh object.
!
lon = 0.0
lat = 0.0
DO i = 1, SIZE(verts)
    CALL global_mesh%get_vert_coords(verts(i), vert_coords)
    lat = lat + vert_coords(2)
    lon = lon + vert_coords(1)
END DO
lat = lat / REAL(SIZE(verts))
lon = lon / REAL(SIZE(verts))

NULLIFY(mesh,global_mesh)

DEALLOCATE(verts)

END SUBROUTINE get_lfric_mesh_coords

END MODULE lfricinp_get_latlon_mod
