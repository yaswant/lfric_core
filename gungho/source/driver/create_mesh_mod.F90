!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
!> @brief Set up and destroy partitioned 3D mesh(es)
!> @details Contains routines to: i) read global ugrid meshes and set up
!>          partitioned 3D mesh(es) ii) destroy partitioned mesh(es)
module create_mesh_mod

  use base_mesh_config_mod,       only: filename, prime_mesh_name, geometry, &
                                        geometry_spherical
  use constants_mod,              only: i_def, str_def, l_def, r_def
  use extrusion_mod,              only: extrusion_type, uniform_extrusion_type
  use extrusion_config_mod,       only: domain_top
  use finite_element_config_mod,  only: cellshape,          &
                                        key_from_cellshape, &
                                        cellshape_triangle, &
                                        cellshape_quadrilateral
  use global_mesh_mod,            only: global_mesh_type
  use global_mesh_collection_mod, only: global_mesh_collection
  use gungho_extrusion_mod,       only: create_extrusion, create_shifted_extrusion
  use create_multigrid_mesh_mod,  only: init_multigrid_mesh
  use log_mod,                    only: log_event,         &
                                        log_scratch_space, &
                                        LOG_LEVEL_INFO,    &
                                        LOG_LEVEL_ERROR
  use mesh_collection_mod,        only: mesh_collection_type, mesh_collection
  use mesh_mod,                   only: mesh_type
  use multigrid_config_mod,       only: l_multigrid
  use ncdf_quad_mod,              only: ncdf_quad_type
  use partition_mod,              only: partition_type,                 &
                                        partitioner_interface,          &
                                        partitioner_cubedsphere_serial, &
                                        partitioner_cubedsphere,        &
                                        partitioner_planar
  use partitioning_config_mod,    only: panel_decomposition,        &
                                        panel_xproc, panel_yproc,   &
                                        PANEL_DECOMPOSITION_AUTO,   &
                                        PANEL_DECOMPOSITION_ROW,    &
                                        PANEL_DECOMPOSITION_COLUMN, &
                                        PANEL_DECOMPOSITION_CUSTOM

  use subgrid_config_mod,         only: dep_pt_stencil_extent, &
                                        rho_approximation_stencil_extent
  use transport_config_mod,       only: scheme, operators,    &
                                        fv_flux_order,        &
                                        fv_advective_order,   &
                                        operators_fv,         &
                                        scheme_horz_cosmic,   &
                                        scheme_yz_bip_cosmic, &
                                        scheme_cosmic_3D
  use mixing_config_mod,          only: smagorinsky

  use ugrid_2d_mod,               only: ugrid_2d_type
  use ugrid_file_mod,             only: ugrid_file_type

  implicit none

  private
  public :: init_mesh, final_mesh

contains

!> @brief Generates a mesh and determines the basis functions and dofmaps
!> @details This will be replaced with code that reads the information in
!> @param[in] local_rank        Number of the MPI rank of this process
!> @param[in] total_ranks       Total number of MPI ranks in this job
!> @param[out] prime_mesh_id    Mesh id of partitioned prime mesh
!> @param[out] twod_mesh_id     Mesh id of the 2D (surface) mesh
!> @param[out] shifted_mesh_id  Mesh id of vertically shifted mesh with an extra level
subroutine init_mesh( local_rank, total_ranks, prime_mesh_id, &
                      twod_mesh_id, shifted_mesh_id )

  implicit none

  integer(i_def), intent(in)  :: local_rank
  integer(i_def), intent(in)  :: total_ranks
  integer(i_def), intent(out) :: prime_mesh_id
  integer(i_def), intent(out) :: twod_mesh_id
  integer(i_def), intent(out), optional :: shifted_mesh_id

  ! Parameters
  integer(i_def), parameter :: max_factor_iters = 10000

  ! Local variables
  procedure (partitioner_interface), pointer :: partitioner_ptr => null()
  type(global_mesh_type),            pointer :: global_mesh_ptr => null()
  class(extrusion_type),         allocatable :: extrusion
  class(extrusion_type),         allocatable :: shifted_extrusion
  type(partition_type)                       :: partition
  type(uniform_extrusion_type)               :: extrusion_sl

  ! max_stencil_depth is the maximum depth (of cells outside the cell over
  ! which the stencil is based) of the stencil to be used on fields with
  ! this partition.
  !
  ! A single cell stencil will, therefore, have a  max_stencil_depth=0.
  ! A nine-point square region stencil will have max_stencil_depth=1
  !
  !> @todo max_stencil_depth will eventually become either a configuration
  !>       item, or will be autogenerated from kernel metadata, but for now
  !>       it is just hard-coded
  integer(i_def) :: max_stencil_depth

  ! Number of ranks the mesh is partitioned over in the x- and y-directions
  ! (across a single face for a cubed-sphere mesh)
  integer(i_def) :: xproc, yproc

  integer(i_def) :: ranks_per_panel
  integer(i_def) :: start_factor
  integer(i_def) :: end_factor
  integer(i_def) :: fact_count
  integer(i_def) :: global_mesh_id
  integer(i_def) :: npanels
  integer(i_def) :: max_fv_stencil
  integer(i_def) :: n_meshes
  integer(i_def) :: i
  logical(l_def) :: found_factors

  character(str_def)  :: domain_desc
  character(str_def)  :: partition_desc
  type(ugrid_2d_type) :: ugrid_2d

  class(ugrid_file_type), allocatable :: file_handler
  character(str_def),     allocatable :: mesh_names(:)


  allocate( mesh_collection, &
            source=mesh_collection_type() )

  call log_event( "Setting up partition mesh(es)", LOG_LEVEL_INFO )

  ! Currently only quad elements are fully functional
  if (cellshape /= cellshape_quadrilateral) then
    call log_event( "Reference_element must be QUAD for now...", &
                    LOG_LEVEL_ERROR )
  end if

  ! Setup the partitioning strategy
  if (geometry == geometry_spherical) then

    npanels = 6
    if (total_ranks == 1 .or. mod(total_ranks,6) == 0) then
      ranks_per_panel = total_ranks/6
      domain_desc = "6x"
    else
      call log_event( "Total number of processors must be a "//   &
                      "multiple of 6 for a cubed-sphere domain.", &
                      LOG_LEVEL_ERROR )
    end if

    if (total_ranks == 1) then
      if (scheme == scheme_horz_cosmic) then
        call log_event( "For Cosmic the total number of processors must be "// &
                        "greater than 1 and a multiple of 6 for a          "// &
                        "cubed-sphere domain.", &
                        LOG_LEVEL_ERROR )
      end if
      ranks_per_panel = 1
      partitioner_ptr => partitioner_cubedsphere_serial
      call log_event( "Using serial cubed sphere partitioner", &
                      LOG_LEVEL_INFO )
    else
      partitioner_ptr => partitioner_cubedsphere
      call log_event( "Using parallel cubed sphere partitioner", &
                      LOG_LEVEL_INFO )
    end if
  else
    npanels = 1
    ranks_per_panel = total_ranks
    domain_desc = ""

    partitioner_ptr => partitioner_planar
    call log_event( "Using planar mesh partitioner ", &
                    LOG_LEVEL_INFO )
  end if

  select case(panel_decomposition)

  case( PANEL_DECOMPOSITION_AUTO )

    ! For automatic partitioning, try to partition into the squarest
    ! possible partitions by finding the two factors of ranks_per_panel
    ! that are closest to sqrt(ranks_per_panel). If two factors can't
    ! be found after max_factor_iters attempts, they would provide
    ! partitions that are too un-square, so an error is produced.
    start_factor  = nint(sqrt(real(ranks_per_panel, kind=r_def)), kind=i_def)
    end_factor    = max(1,(start_factor-max_factor_iters))
    found_factors = .false.
    do fact_count = start_factor, end_factor, -1
      if (mod(ranks_per_panel,fact_count) == 0) then
        found_factors = .true.
        exit
      end if
    end do

    if (found_factors) then
      xproc = fact_count
      yproc = ranks_per_panel/fact_count
    else
      call log_event( "Could not automatically partition domain.", &
                      LOG_LEVEL_ERROR )
    end if

  case( PANEL_DECOMPOSITION_ROW )
    xproc = ranks_per_panel
    yproc = 1

  case( PANEL_DECOMPOSITION_COLUMN )
    xproc = 1
    yproc = ranks_per_panel

  case( PANEL_DECOMPOSITION_CUSTOM )
    ! Use the values provided from the partitioning namelist
    xproc = panel_xproc
    yproc = panel_yproc

    if (xproc*yproc /= ranks_per_panel) then
      call log_event( "The values of panel_xproc and panel_yproc "// &
                      "are inconsistent with the total number of "// &
                      "processors available.", LOG_LEVEL_ERROR )
    end if

  case default

    call log_event( "Missing entry for panel decomposition, "// &
                    "specify 'auto' if unsure.", LOG_LEVEL_ERROR )

  end select


  write(log_scratch_space, '(2(A,I0),A)' )   &
        'Using ', xproc*yproc*npanels, ' '// &
        'partition(s) across ', npanels, ' domain'
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  if (total_ranks > 1) then
    write(log_scratch_space, '(2(A,I0))' )   &
        'Panel decomposition: ', xproc,'x', yproc
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  end if

  ! Determine max_stencil_depth
  max_stencil_depth = 1
  ! Smagorinsky (or boundary layers) appears to need larger haloes
  if ( smagorinsky ) max_stencil_depth = max(max_stencil_depth,2)
  if (operators == operators_fv) then
    ! Need larger haloes for fv operators
    max_fv_stencil = max(fv_flux_order,fv_advective_order)/2
    max_stencil_depth = max(max_stencil_depth,max_fv_stencil)
  end if

  if (scheme == scheme_yz_bip_cosmic .or. &
      scheme == scheme_horz_cosmic   .or. &
      scheme == scheme_cosmic_3D   ) then
    max_stencil_depth = max(max_stencil_depth,      &
                            dep_pt_stencil_extent + &
                            rho_approximation_stencil_extent)
  end if

  ! Interrogate ugrid file to get the names of all the
  ! contained mesh topologies
  allocate( ncdf_quad_type :: file_handler )
  call ugrid_2d%set_file_handler( file_handler )
  call ugrid_2d%get_n_meshes( trim(filename), n_meshes )

  allocate( mesh_names(n_meshes) )
  call ugrid_2d%get_mesh_names( trim(filename), mesh_names )

  ! Read in prime mesh only
  ! Other meshes may need to be read in when multiple meshes
  ! are use i.e. in function space chains.
  do i=1, n_meshes
    if (trim(mesh_names(i)) == trim(prime_mesh_name)) then
      global_mesh_id = global_mesh_collection %                 &
                           add_new_global_mesh ( filename,      &
                                                 mesh_names(i), &
                                                 npanels )

      global_mesh_ptr => global_mesh_collection % &
                             get_global_mesh( global_mesh_id )
      exit
    end if
  end do

  if (.not. associated(global_mesh_ptr)) &
      call log_event( "Global mesh not in collection", LOG_LEVEL_ERROR )

  ! Generate the partition object
  partition = partition_type( global_mesh_ptr,   &
                              partitioner_ptr,   &
                              xproc,             &
                              yproc,             &
                              max_stencil_depth, &
                              local_rank,        &
                              total_ranks )

  allocate(extrusion, source=create_extrusion() )

  ! Generate the mesh
  call log_event( "Creating prime mesh", LOG_LEVEL_INFO )
  prime_mesh_id = mesh_collection%add_new_mesh( global_mesh_ptr,  &
                                                partition,        &
                                                extrusion )

  write(log_scratch_space,'(A,I0,A)') &
      "Prime mesh created (id:", prime_mesh_id, ")"
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  ! Generate a '2d' mesh
  ! probably only works for cartesian domains, as atmos_bottom hard-wired
  ! to 0 currently...
  extrusion_sl = uniform_extrusion_type( 0.0_r_def, 1.0_r_def, 1_i_def )
  twod_mesh_id = mesh_collection%add_new_mesh( global_mesh_ptr,         &
                                               partition,               &
                                               extrusion_sl )

  if (l_multigrid) then
    ! Done this way as e do not currently have multiple global
    ! meshes and their associated maps in a single ugrid file
    call init_multigrid_mesh( prime_mesh_id,           &
                              twod_mesh_id,            &
                              partitioner_ptr,         &
                              xproc, yproc,            &
                              max_stencil_depth,       &
                              local_rank, total_ranks, &
                              npanels )
  end if


  if (present(shifted_mesh_id)) then
    allocate(shifted_extrusion, source=create_shifted_extrusion(extrusion) )

    call log_event( "Creating shifted mesh", LOG_LEVEL_INFO )

    shifted_mesh_id = mesh_collection%add_new_mesh( global_mesh_ptr,  &
                                                    partition,        &
                                                    shifted_extrusion )

    deallocate(shifted_extrusion)
  end if

  deallocate(extrusion)

  return
end subroutine init_mesh

!> @brief Finalises the mesh_collection
subroutine final_mesh()

  implicit none

  if (allocated(mesh_collection)) then
    call mesh_collection%clear()
    deallocate(mesh_collection)
  end if

end subroutine final_mesh

end module create_mesh_mod
