!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Defines various mesh details.

!> @details information about the mesh is held in  here.
module mesh_mod
use constants_mod,  only : r_def
implicit none

!> Total number of horizontal cells in the domain on the local partition 
integer :: num_cells
!> For a biperiodic mesh, the number of horizontal cells in the x-direction of the global domain
!> For a cubed-sphere mesh it is the number of cells across a face.
integer :: num_cells_x
!> Number of horizontal cells in the y-direction of the global domain (not used for cubed-sphere meshes)
integer :: num_cells_y
!> Number of vertical layers
integer :: num_layers
!> Order of the function space
integer :: element_order
!> Flag for whether mesh is on a sphere or not
logical :: l_spherical

!> Number of unique dofs in a particular function space (4,:), either globally (:,1) or per cell (:,2)
integer :: w_unique_dofs(4,2)
!> Number of dofs in a particular function space (4,:) per entity (:,0:3)
integer :: w_dof_entity(4,0:3)

!> Grid spacing in the x-direction
real(kind=r_def)  :: dx
!> Grid spacing in the y-direction
real(kind=r_def)  :: dy
!> Grid spacing in the z-direction
real(kind=r_def)  :: dz

!> Total number of MPI ranks available
integer :: total_ranks
!> The local rank number
integer :: local_rank

!> Number of processors along x-direction
integer :: xproc
!> Number of processors along y-direction
integer :: yproc

!> Global ids of the cells on this partition
integer, allocatable :: partitioned_cells( : )

!> Number of cells that are wholly owned by the partition (i.e. all dofs in these cells are wholly owned by the partition)
integer :: num_core
!> Number of cells that are owned by the partition, but may have dofs that are also owned by halo cells
integer :: num_owned
!> Number of cells in the halo for this partition
integer :: num_halo

end module mesh_mod

