!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Sets up data required for the the function spaces

!> @details This code generates a mesh and determines the basis functions and
!> dofmaps. This will be replaced with code that reads this in from a mesh
!> generation and paritioning pre-processor stage.

! There are no tests for this code as this will be replaced.

module set_up_mod

  use base_mesh_config_mod,      only : filename, &
                                        geometry, &
                                        base_mesh_geometry_spherical
  use constants_mod,             only : i_def, r_def, str_def, PI
  use extrusion_config_mod,      only : number_of_layers,           &
                                        domain_top,                 &
                                        method,                     &
                                        extrusion_method_uniform,   &
                                        extrusion_method_quadratic, &
                                        extrusion_method_geometric, &
                                        extrusion_method_dcmip
  use finite_element_config_mod, only : finite_element_shape_triangle, &
                                        finite_element_shape_quadrilateral
  use global_mesh_mod,           only : global_mesh_type
  use reference_element_mod,     only : reference_cube, reference_element, &
                                        nfaces, nedges, nverts
  use mesh_mod,                  only : mesh_type
  use log_mod,                   only : log_event, LOG_LEVEL_INFO
  use partition_mod,             only : partition_type,                   &
                                        partitioner_interface,            &
                                        partitioner_cubedsphere_serial,   &
                                        partitioner_cubedsphere,          &
                                        partitioner_biperiodic

  implicit none

contains

!> @brief Generates a mesh and determines the basis functions and dofmaps
!> @details This will be replaced with code that reads the information in
!> @param[out] mesh Mesh object to run model on
!> @param[in] local_rank Number of the MPI rank of this process
!> @param[in] total_ranks Total number of MPI ranks in this job
  subroutine set_up(mesh, local_rank, total_ranks)

    implicit none

    type (mesh_type), intent(out), target   :: mesh
    integer(i_def), intent(in) :: local_rank
    integer(i_def), intent(in) :: total_ranks

    type (global_mesh_type)    :: global_mesh
    type (partition_type)      :: partition


    procedure (partitioner_interface), pointer :: partitioner_ptr => null ()

    ! Number of ranks the mesh is partitioned over in the x- and y-directions
    ! (across a single face for a cubed-sphere mesh)
    integer(i_def) :: xproc, yproc


    ! Until we have configuration, set up xproc and yproc for square decomposition
    xproc=max(1,int(sqrt(real(total_ranks/6.0))))
    yproc=max(1,int(sqrt(real(total_ranks/6.0))))
!> @todo Eventually xproc and yproc will be inputted into Dynamo (and not hard-coded).
!>       When this happens their values will need to be checked to make sure they are
!>       sensible



    call log_event( "set_up: Generating/reading the mesh", LOG_LEVEL_INFO )

    ! Currently only quad elements are fully functional
    if ( reference_element /= finite_element_shape_quadrilateral ) then
      call log_event( "set_up: Reference_element must be QUAD for now...", LOG_LEVEL_INFO )
    end if
    ! Setup reference cube
    call reference_cube()

    ! Generate the global mesh and choose a partitioning strategy by setting
    ! a function pointer to point at the appropriate partitioning routine
    global_mesh = global_mesh_type( filename )
    if ( geometry == base_mesh_geometry_spherical ) then
      partitioner_ptr => partitioner_cubedsphere_serial
      call log_event( "set_up: Setting up cubed sphere partitioner ", LOG_LEVEL_INFO )
    else
      partitioner_ptr => partitioner_biperiodic
      call log_event( "set_up: Setting up biperiodic plane partitioner ", LOG_LEVEL_INFO )
    end if

    ! Generate the partition object
    partition = partition_type( global_mesh,    &
                               partitioner_ptr, &
                               xproc,           &
                               yproc,           &
                               1,               &
                               local_rank,      &
                               total_ranks)



    ! Generate the mesh
    mesh =  mesh_type ( partition, global_mesh, &
                        number_of_layers, domain_top, method )
    call mesh%set_colours()


    return
  end subroutine set_up
end module set_up_mod
