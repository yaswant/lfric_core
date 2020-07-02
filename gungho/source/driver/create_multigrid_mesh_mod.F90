!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module create_multigrid_mesh_mod

use base_mesh_config_mod,       only: prime_mesh_name
use constants_mod,              only: i_def, r_def, str_def
use extrusion_mod,              only: extrusion_type, uniform_extrusion_type
use global_mesh_mod,            only: global_mesh_type
use global_mesh_collection_mod, only: global_mesh_collection
use gungho_extrusion_mod,       only: create_extrusion
use log_mod,                    only: log_event, log_scratch_space,    &
                                      LOG_LEVEL_INFO, LOG_LEVEL_TRACE, &
                                      LOG_LEVEL_ERROR
use mesh_collection_mod,        only: mesh_collection
use mesh_mod,                   only: mesh_type
use multigrid_config_mod,       only: ugrid, multigrid_chain_nitems
use partition_mod,              only: partition_type, partitioner_interface

implicit none

private
public :: init_multigrid_mesh, mesh_ids, twod_mesh_ids

integer(i_def), allocatable :: mesh_ids(:), twod_mesh_ids(:)

contains

subroutine init_multigrid_mesh( prime_mesh_id,           &
                                twod_mesh_id,            &
                                partitioner,             &
                                xproc, yproc,            &
                                max_stencil_depth,       &
                                local_rank, total_ranks, &
                                npanels                  )
implicit none


integer, intent(in) :: prime_mesh_id
integer, intent(in) :: twod_mesh_id
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
class(extrusion_type), allocatable :: extrusion
type(uniform_extrusion_type)       :: extrusion_sl

character(str_def) :: mesh_name

mesh_name = prime_mesh_name
global_mesh_ids(:) = 0
allocate(mesh_ids(multigrid_chain_nitems))
allocate(twod_mesh_ids(multigrid_chain_nitems))

mesh_ids(1) = prime_mesh_id
prime_mesh  => mesh_collection%get_mesh(prime_mesh_id)
global_mesh_ids(1) = prime_mesh%get_global_mesh_id()
call global_mesh_collection % set_next_source_mesh(global_mesh_ids(1))

allocate( extrusion, source=create_extrusion() )
extrusion_sl = uniform_extrusion_type( 0.0_r_def, 1.0_r_def, 1_i_def )

if (multigrid_chain_nitems > size(ugrid, 1)) then
  write(log_scratch_space, &
        '(A, I0, A, I0)') 'Requested ',                                       &
                          multigrid_chain_nitems,                             &
                          ' multigrid levels but only specified meshes for ', &
                          size(ugrid, 1)
  call log_event( log_scratch_space, log_level_error )
end if

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

  mesh_ids(i) = mesh_collection % add_new_mesh( global_mesh, &
                                                partition,   &
                                                extrusion )

end do

twod_mesh_ids(1) = twod_mesh_id
prime_mesh => mesh_collection%get_mesh( twod_mesh_id )
do i=2, multigrid_chain_nitems

  global_mesh => global_mesh_collection % get_global_mesh( global_mesh_ids(i) )

  partition =  partition_type( global_mesh, partitioner,        &
                               xproc, yproc, max_stencil_depth, &
                               local_rank, total_ranks )

  twod_mesh_ids(i) = mesh_collection % add_new_mesh( global_mesh, &
                                                     partition,   &
                                                     extrusion_sl )

end do

deallocate(extrusion)

return
end subroutine init_multigrid_mesh

end module create_multigrid_mesh_mod
