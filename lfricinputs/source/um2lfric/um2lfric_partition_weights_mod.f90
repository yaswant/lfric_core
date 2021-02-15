! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE um2lfric_partition_weights_mod

! LFRicInputs modules
USE lfricinp_regrid_weights_type_mod, ONLY: lfricinp_regrid_weights_type
! LFRic modules
USE field_mod,       ONLY: lfric_field_type => field_type

IMPLICIT NONE

PRIVATE

PUBLIC :: um2lfric_partition_weights

CONTAINS

SUBROUTINE um2lfric_set_local_weights(weights, field)

! Description: Set local weights indices so that the weights type on each partition
!  only contains indices relevant to that partition.

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int32, real64

! LFRic modules
USE mesh_mod,        ONLY: mesh_type
USE local_mesh_mod,  ONLY: local_mesh_type
IMPLICIT NONE

! Arguments

TYPE(lfricinp_regrid_weights_type), INTENT(INOUT) :: weights
TYPE(lfric_field_type), INTENT(IN)       :: field

! Local variables
TYPE(mesh_type), POINTER :: mesh
TYPE(local_mesh_type), POINTER :: local_mesh

! Local partition versions of weights arrays
INTEGER(KIND=int32), ALLOCATABLE :: dst_address_local_tmp(:)
INTEGER(KIND=int32), ALLOCATABLE :: src_address_local_tmp(:)
INTEGER(KIND=int32), ALLOCATABLE :: src_address_local_2D_tmp(:,:)
REAL(KIND=real64), ALLOCATABLE   :: remap_matrix_local_tmp(:,:)

INTEGER(KIND=int32) :: local_links, i_link, w

mesh => field % get_mesh()
local_mesh => mesh % get_local_mesh()

! Initially allocate to be the same size as global arrays
ALLOCATE( dst_address_local_tmp( weights%num_links ) )
ALLOCATE( src_address_local_tmp( weights%num_links ) )
ALLOCATE( src_address_local_2D_tmp( weights%num_links, 2) )
ALLOCATE( remap_matrix_local_tmp( weights%num_wgts, weights%num_links) )
! Convert weights from global indices to the local partitions id/indices
! Use local_links to count how many links on local partition
local_links = 0
DO i_link = 1, weights%num_links
  ! The call to get_lid_from_gid will return -1 if the global id is not
  ! known to this partition. Also the partition will also contain halo cells.
  ! However, we only want to select cells wholy owned by this partition. So reject
  ! any addresses that are -1 or greater than ncells_2D
  IF ( local_mesh % get_lid_from_gid(weights%dst_address(i_link)) /= -1 .AND.  &
       local_mesh % get_lid_from_gid(weights%dst_address(i_link)) <=          &
         mesh % get_ncells_2d() ) THEN
    local_links = local_links + 1
    dst_address_local_tmp(local_links) = &
                      local_mesh % get_lid_from_gid(weights%dst_address(i_link))
    ! Weights and source address also need filtering
    src_address_local_tmp(local_links) = weights%src_address(i_link)
    src_address_local_2D_tmp(local_links, 1) = weights%src_address_2D(i_link, 1)
    src_address_local_2D_tmp(local_links, 2) = weights%src_address_2D(i_link, 2)
    DO w = 1, weights%num_wgts
      remap_matrix_local_tmp(w, local_links) = weights%remap_matrix(w, i_link)
    END DO
  END IF
END DO

! Deallocate and reallocate the arrays in the main weights type using
! the local sizes as calculated by the local_links counter in the above loop
IF (ALLOCATED(weights%src_address_2D)) DEALLOCATE(weights%src_address_2D)
ALLOCATE(weights%src_address_2D(local_links, 2))
IF (ALLOCATED(weights%src_address)) DEALLOCATE(weights%src_address)
ALLOCATE( weights%src_address( local_links ))
IF (ALLOCATED(weights%dst_address)) DEALLOCATE(weights%dst_address)
ALLOCATE( weights%dst_address( local_links ))
IF (ALLOCATED(weights%remap_matrix)) DEALLOCATE(weights%remap_matrix)
ALLOCATE(weights%remap_matrix( weights%num_wgts, local_links))

weights % num_links = local_links
! Copy values from the tmp versions of arrays to the main weights type
DO i_link = 1, weights%num_links
  weights%dst_address(i_link) = dst_address_local_tmp(i_link)
  weights%src_address(i_link) = src_address_local_tmp(i_link)
  weights%src_address_2D(i_link,1) = src_address_local_2D_tmp(i_link,1)
  weights%src_address_2D(i_link,2) = src_address_local_2D_tmp(i_link,2)
  weights%remap_matrix(:,i_link) = remap_matrix_local_tmp(:,i_link)
END DO

END SUBROUTINE um2lfric_set_local_weights

!-------------------------------------------------------

SUBROUTINE um2lfric_partition_weights()
USE lfricinp_stashmaster_mod, ONLY: stashcode_theta, stashcode_u, &
    stashcode_v
USE um2lfric_regrid_weights_mod, ONLY: get_weights
USE lfricinp_lfric_driver_mod, ONLY: lfric_fields
USE lfricinp_stash_to_lfric_map_mod, ONLY: get_field_name

IMPLICIT NONE

! Description: Control routine to call weights partitioning for each grid type
! Will need additional logic once we have a dynamic list of fields
TYPE(lfricinp_regrid_weights_type), POINTER :: weights
TYPE(lfric_field_type), POINTER    :: field

! Get an LFRic field that matches the function space of the weights files.
! Currently we only support cell centres, use theta for now.

field => lfric_fields%get_field(get_field_name(stashcode_theta))

! U Points
weights => get_weights(stashcode_u)
CALL um2lfric_set_local_weights(weights, field)

! V Points
weights => get_weights(stashcode_v)
CALL um2lfric_set_local_weights(weights, field)

! Theta/P points
weights => get_weights(stashcode_theta)
CALL um2lfric_set_local_weights(weights, field)

END SUBROUTINE um2lfric_partition_weights

END MODULE um2lfric_partition_weights_mod
