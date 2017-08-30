!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module init_multigrid_mesh_mod

use constants_mod,              only: i_def, str_def
use global_mesh_collection_mod, only: global_mesh_collection
use mesh_collection_mod,        only: mesh_collection
use mesh_mod,                   only: mesh_type
use global_mesh_mod,            only: global_mesh_type
use partition_mod,              only: partition_type, partitioner_interface
use extrusion_config_mod,       only: number_of_layers, domain_top, method
use base_mesh_config_mod,       only: prime_mesh_name
use log_mod,                    only: log_event, log_scratch_space,    &
                                      LOG_LEVEL_INFO, LOG_LEVEL_TRACE, &
                                      LOG_LEVEL_ERROR

use multigrid_config_mod, only: l_multigrid, ugrid, order, &
                                continuity, multigrid_chain_nitems

implicit none

private
public :: init_multigrid_mesh, mesh_ids

integer(i_def), allocatable :: mesh_ids(:)

contains

subroutine init_multigrid_mesh( prime_mesh_id,           &
                                partitioner,             &
                                xproc, yproc,            &
                                max_stencil_depth,       &
                                local_rank, total_ranks, &
                                npanels                  )
implicit none


integer, intent(in) :: prime_mesh_id
procedure(partitioner_interface),  intent(in), pointer :: partitioner

integer(i_def), intent(in) :: xproc, yproc
integer(i_def), intent(in) :: max_stencil_depth
integer(i_def), intent(in) :: local_rank, total_ranks
integer(i_def), intent(in) :: npanels

integer(i_def) :: global_mesh_ids (multigrid_chain_nitems)
integer(i_def) :: i

type(mesh_type),        pointer :: prime_mesh => null()
type(global_mesh_type), pointer :: global_mesh => null()
type(partition_type) :: partition

character(str_def) :: mesh_name

mesh_name = prime_mesh_name
global_mesh_ids(:) = 0
allocate(mesh_ids(multigrid_chain_nitems))

!----------------------------------------------------------------------------
! First global mesh should be the prime mesh so check that first
if ( trim(ugrid(1)) /= 'prime' ) then
  ! Error illegal input for multigrid
  write( log_scratch_space,'(A)')                                &
      'First function space in function space chains must be ' //&
      'based on the "prime" mesh.'//trim(ugrid(1))
  call log_event( log_scratch_space, LOG_LEVEL_ERROR )
end if


mesh_ids(1) = prime_mesh_id
prime_mesh  => mesh_collection%get_mesh(prime_mesh_id)
global_mesh_ids(1) = prime_mesh%get_global_mesh_id()
call global_mesh_collection % set_next_source_mesh(global_mesh_ids(1))


! Read global meshes into the global mesh collection,
! maps should be created in the order they are added.
!===================================================================
do i=2, multigrid_chain_nitems
  global_mesh_ids(i) = global_mesh_collection %                    &
                                   add_new_global_mesh( ugrid(i),  & 
                                                        mesh_name, &
                                                        npanels )

  global_mesh => global_mesh_collection % get_global_mesh( global_mesh_ids(i) )

  partition =  partition_type( global_mesh, partitioner,        &
                               xproc, yproc, max_stencil_depth, &
                               local_rank, total_ranks )

  mesh_ids(i) = mesh_collection % add_new_mesh               &
                             ( global_mesh, partition,       &
                               number_of_layers, domain_top, &
                               method )

end do

return
end subroutine init_multigrid_mesh

end module init_multigrid_mesh_mod
