!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Controls grid related information used by the model

module gravity_wave_grid_mod

  use base_mesh_config_mod,           only : prime_mesh_name
  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type
  use global_mesh_collection_mod,     only : global_mesh_collection, &
                                             global_mesh_collection_type
  use create_mesh_mod,                only : init_mesh
  use create_fem_mod,                 only : init_fem
  use fs_continuity_mod,              only : W3
  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection_type, &
                                             function_space_collection
  use function_space_chain_mod,       only : function_space_chain_type
  use log_mod,                        only : log_event,          &
                                             log_scratch_space,  &
                                             LOG_LEVEL_INFO
  use formulation_config_mod,         only : l_multigrid
  use mpi_mod,                        only : get_comm_size, get_comm_rank

  implicit none

  private
  public initialise_grid

contains

  !> @brief Initialises grid related information used by the model
  !> @param [inout] mesh_id The identifier of the primary mesh
  !> @param [inout] twod_mesh_id The identifier of the primary 2d mesh
  !> @param [inout] chi Function space representatio of the 3d mesh
  !> @param [inout] multigrid_function_space_chain An ordered  list of function
  !>                                               spaces used by multigrid
  subroutine initialise_grid(mesh_id, twod_mesh_id, chi, &
                             multigrid_function_space_chain)

    implicit none

    integer(i_def),   intent(inout) :: mesh_id, twod_mesh_id
    type(field_type), intent(inout) :: chi(3)
    type(function_space_chain_type), intent(inout) :: &
                                       multigrid_function_space_chain

    type(function_space_type), pointer :: function_space

    integer(i_def)     :: total_ranks, local_rank
    integer(i_def)     :: i

    integer(i_def), allocatable :: multigrid_mesh_ids(:)

    allocate( global_mesh_collection, &
              source = global_mesh_collection_type() )

    ! Get the rank information from the virtual machine
    total_ranks = get_comm_size()
    local_rank  = get_comm_rank()

    ! Create the mesh
    call init_mesh( local_rank, total_ranks, mesh_id,     &
                    twod_mesh_id=twod_mesh_id,            &
                    multigrid_mesh_ids=multigrid_mesh_ids )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh_id, chi )

    if (l_multigrid) then

      do i=1, size(multigrid_mesh_ids)

        ! Make sure this function_space is in the collection
        function_space => &
            function_space_collection%get_fs( multigrid_mesh_ids(i), 0, W3 )

        write( log_scratch_space,"(A,I0,A)")                       &
             'Adding function_space id ', function_space%get_id(), &
             ' to multigrid function_space chain'
        call log_event( log_scratch_space, LOG_LEVEL_INFO )

        call multigrid_function_space_chain%add( function_space )

      end do
    end if

    nullify( function_space )

    ! Full global meshes no longer required, so reclaim
    ! the memory from global_mesh_collection
    if (allocated(global_mesh_collection)) deallocate(global_mesh_collection)

  end subroutine initialise_grid

end module gravity_wave_grid_mod
