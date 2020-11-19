!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
!> @brief Set up and destroy partitioned 3D mesh(es)
!> @details Contains routines to:
!>           i) Read global ugrid meshes and set up partitioned 3D mesh(es)
!>          ii) Destroy partitioned mesh(es)
module create_mesh_mod


  use constants_mod,              only: i_def, str_def, l_def, &
                                        r_def, imdi, i_native
  use extrusion_mod,              only: extrusion_type, &
                                        uniform_extrusion_type
  use global_mesh_mod,            only: global_mesh_type
  use log_mod,                    only: log_event,         &
                                        log_scratch_space, &
                                        LOG_LEVEL_INFO,    &
                                        LOG_LEVEL_ERROR
  use mesh_collection_mod,        only: mesh_collection_type
  use mesh_mod,                   only: mesh_type
  use ncdf_quad_mod,              only: ncdf_quad_type
  use partition_mod,              only: partition_type, &
                                        partitioner_interface
  use ugrid_2d_mod,               only: ugrid_2d_type
  use ugrid_mesh_data_mod,        only: ugrid_mesh_data_type
  use ugrid_file_mod,             only: ugrid_file_type

  implicit none

  private
  public  :: init_mesh, final_mesh

  private :: read_global_meshes,           &
             assign_global_intergrid_maps, &
             set_partition_parameters,     &
             create_3d_mesh_partitions,    &
             create_3d_mesh_partition,     &
             add_intergrid_maps,           &
             load_mesh

contains


!===============================================================================
!> @brief Generates a mesh and determines the basis functions and dofmaps
!> @details This will be replaced with code that reads the information in
!>
!> @param[in]  local_rank   Number of the MPI rank of this process
!> @param[in]  total_ranks  Total number of MPI ranks in this job
!===============================================================================
!> @brief Generates a mesh and determines the basis functions and dofmaps
!> @details This will be replaced with code that reads the information in
!> @param[in]  local_rank            Number of the MPI rank of this process
!> @param[in]  total_ranks           Total number of MPI ranks in this job
!> @param[out] prime_mesh_id         Mesh id of partitioned prime mesh
!> @param[out] twod_mesh_id          Mesh id of the 2D (surface) mesh
!> @param[out] shifted_mesh_id       Mesh id of vertically shifted mesh
!>                                   with an extra level
!> @param[out] double_level_mesh_id  Mesh id of vertically double level mesh
!> @param[out] multigrid_mesh_ids    Gungho multigrid chain mesh ids
!> @param[out] multigrid_2D_mesh_ids Gungho multigrid chain 2D_mesh ids

subroutine init_mesh( local_rank, total_ranks, &
                      mesh_id,                 &
                      twod_mesh_id,            &
                      shifted_mesh_id,         &
                      double_level_mesh_id,    &
                      multigrid_mesh_ids,      &
                      multigrid_2D_mesh_ids )

  use base_mesh_config_mod,       only: prime_mesh_name
  use finite_element_config_mod,  only: cellshape,          &
                                        key_from_cellshape, &
                                        cellshape_triangle, &
                                        cellshape_quadrilateral
  use formulation_config_mod,     only: l_multigrid
  use mesh_collection_mod,        only: mesh_collection
  use multigrid_config_mod,       only: chain_mesh_tags

  implicit none

  integer(i_def), intent(in)  :: local_rank
  integer(i_def), intent(in)  :: total_ranks
  integer(i_def), intent(out) :: mesh_id

  integer(i_def), intent(out), optional :: twod_mesh_id
  integer(i_def), intent(out), optional :: shifted_mesh_id
  integer(i_def), intent(out), optional :: double_level_mesh_id
  integer(i_def), intent(out), optional, allocatable :: multigrid_mesh_ids(:)
  integer(i_def), intent(out), optional, allocatable :: multigrid_2d_mesh_ids(:)

  ! Parameters
  integer(i_def), parameter :: max_factor_iters = 10000

  ! Local variables
  integer(i_def) :: xproc  ! Processor ranks in mesh panel x-direction
  integer(i_def) :: yproc  ! Processor ranks in mesh panel y-direction

  ! max_stencil_depth is the maximum depth (of cells outside the cell over
  ! which the stencil is based) of the stencil to be used on fields with
  ! this partition.
  !
  ! A single cell stencil will, therefore, have a max_stencil_depth=0.
  ! A nine-point square region stencil will have max_stencil_depth=1
  !
  !> @todo max_stencil_depth will eventually become either a configuration
  !>       item, or will be autogenerated from kernel metadata, but for now
  !>       it is just hard-coded
  integer(i_def) :: max_stencil_depth

  procedure(partitioner_interface), pointer :: partitioner_ptr => null()

  logical(l_def) :: create_2d_mesh             = .false.
  logical(l_def) :: create_shifted_mesh        = .false.
  logical(l_def) :: create_double_level_mesh   = .false.
  logical(l_def) :: create_multigrid_meshes    = .false.
  logical(l_def) :: create_multigrid_2d_meshes = .false.

  integer(i_def) :: i, n_chain_meshes

  character(str_def) :: mesh_name
  character(str_def) :: mesh_name_A
  character(str_def) :: mesh_name_B

  allocate( mesh_collection, &
            source=mesh_collection_type() )

  call log_event( "Setting up partition mesh(es)", LOG_LEVEL_INFO )

  ! Currently only quad elements are fully functional
  if (cellshape /= cellshape_quadrilateral) then
    call log_event( "Reference_element must be QUAD for now...", &
                    LOG_LEVEL_ERROR )
  end if

  ! Set control flags
  !=================================================================
  if ( present(twod_mesh_id) )         create_2d_mesh           = .true.
  if ( present(shifted_mesh_id) )      create_shifted_mesh      = .true.
  if ( present(double_level_mesh_id) ) create_double_level_mesh = .true.

  if ( l_multigrid ) then
    if ( present(multigrid_mesh_ids) ) then
      create_multigrid_meshes    = .true.
      if (allocated(multigrid_mesh_ids)) deallocate(multigrid_mesh_ids)
      allocate(multigrid_mesh_ids(size(chain_mesh_tags)))
    end if
    if ( present(multigrid_2d_mesh_ids) )then
      create_multigrid_2d_meshes = .true.
      if (allocated(multigrid_2d_mesh_ids)) deallocate(multigrid_2d_mesh_ids)
      allocate(multigrid_2d_mesh_ids(size(chain_mesh_tags)))
    end if
  end if


  ! 1.0 Read in all the required global meshes from the mesh file.
  !=================================================================
  call read_global_meshes()


  ! 2.0 Read in and assign any intergrid maps connected with
  !     global meshes in the global mesh collection.
  !=================================================================
  call assign_global_intergrid_maps()


  ! 3.0 Set runtime constant partition parameters.
  !=================================================================
  call set_partition_parameters( total_ranks,       &
                                 xproc, yproc,      &
                                 max_stencil_depth, &
                                 partitioner_ptr )

  ! 4.0 Partition global meshes and extrude into 3D local meshes
  !     after reading in ALL the global meshes required for the
  !     run configuration.
  !=================================================================
  call create_3D_mesh_partitions( local_rank, total_ranks,  &
                                  xproc, yproc,             &
                                  max_stencil_depth,        &
                                  partitioner_ptr,          &
                                  create_2d_mesh,           &
                                  create_shifted_mesh,      &
                                  create_double_level_mesh, &
                                  create_multigrid_meshes )

  ! 5.0 Assign maps to local 3D meshes
  !=================================================================
  n_chain_meshes = size(chain_mesh_tags)
  if (create_multigrid_meshes) then
    do i=1, n_chain_meshes-1
      mesh_name_A = chain_mesh_tags(i)
      mesh_name_B = chain_mesh_tags(i+1)
      call add_intergrid_maps( mesh_name_A, mesh_name_B )
      multigrid_mesh_ids(i) = mesh_collection%get_mesh_id(mesh_name_A)
    end do
    multigrid_mesh_ids(n_chain_meshes) = mesh_collection%get_mesh_id(mesh_name_B)
  end if

  if (create_multigrid_2d_meshes) then
    do i=1, n_chain_meshes-1
      mesh_name_A = trim(chain_mesh_tags(i))//'_2d'
      mesh_name_B = trim(chain_mesh_tags(i+1))//'_2d'
      call add_intergrid_maps( mesh_name_A, mesh_name_B )
      multigrid_2d_mesh_ids(i) = mesh_collection%get_mesh_id(mesh_name_A)
    end do
    multigrid_2d_mesh_ids(n_chain_meshes) = mesh_collection%get_mesh_id(mesh_name_B)
  end if


  ! 6.0 Extract out mesh ids
  !=================================================================
  mesh_name = prime_mesh_name
  mesh_id = mesh_collection%get_mesh_id(mesh_name)

  if ( create_2d_mesh ) then
    mesh_name = trim(prime_mesh_name)//'_2d'
    twod_mesh_id = mesh_collection%get_mesh_id(mesh_name)
  end if

  if ( create_shifted_mesh ) then
    mesh_name = trim(prime_mesh_name)//'_shifted'
    shifted_mesh_id = mesh_collection%get_mesh_id(mesh_name)
  end if

  if ( create_double_level_mesh ) then
    mesh_name = trim(prime_mesh_name)//'_double'
    double_level_mesh_id = mesh_collection%get_mesh_id(mesh_name)
  end if

  return
end subroutine init_mesh


!>==============================================================================
!> @brief   Reads in global meshes from file.
!>          (private subroutine)
!> @details All global meshes are loaded in similar fashion. The meshes
!>          to be loaded are determined by the model configuration.
!>==============================================================================
subroutine read_global_meshes()

  use base_mesh_config_mod,   only: prime_mesh_name, &
                                    geometry, geometry_spherical
  use formulation_config_mod, only: l_multigrid
  use multigrid_config_mod,   only: chain_mesh_tags

  implicit none

  integer(i_def) :: n_panels

  if (geometry == geometry_spherical) then
    n_panels = 6
  else
    n_panels = 1
  end if

  write(log_scratch_space, '(A,I0,A)' )        &
      'Creating global meshes comprising of ', &
      n_panels, ' domain(s)'
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  ! 1.0 Read in prime mesh first by default
  !----------------------------------------------
  write(log_scratch_space,'(A)') &
      'Reading prime global mesh: "'//trim(prime_mesh_name)//'"'
  call log_event(log_scratch_space, LOG_LEVEL_INFO)

  call load_mesh([prime_mesh_name], n_panels)

  ! 2.0 Read in any other global meshes required
  !     by other additional configuraed schemes
  !----------------------------------------------
  if (l_multigrid) call load_mesh(chain_mesh_tags, n_panels)

  return
end subroutine read_global_meshes


!>==============================================================================
!> @brief   Reads in and assigns available global intergrid maps from file.
!>          (private subroutine)
!> @details Global meshes which have been read into the models global mesh
!>          collection will have a list of target mesh names. These target mesh
!>          names (if any) indicate the valid intergrid maps avaiable in the
!>          mesh file.
!>
!>          This routine will read in in the appropriate intergrid maps and
!>          assign them to the correct global mesh object.
!>==============================================================================
subroutine assign_global_intergrid_maps()

  use base_mesh_config_mod,       only: filename
  use global_mesh_collection_mod, only: global_mesh_collection

  implicit none

  type(ncdf_quad_type) :: file_handler

  character(str_def), allocatable :: global_mesh_names(:)
  character(str_def), allocatable :: target_mesh_names(:)
  integer(i_def),     allocatable :: mesh_map(:,:)

  integer(i_def) :: i, j
  integer(i_def) :: n_meshes

  type(global_mesh_type), pointer :: global_mesh =>null()
  type(global_mesh_type), pointer :: target_mesh =>null()

  ! Read in the maps for each global mesh
  !=================================================================
  call file_handler%file_open(trim(filename))

  global_mesh_names = global_mesh_collection%get_mesh_names()
  n_meshes = global_mesh_collection%n_meshes()

  do i=1, n_meshes
    global_mesh => global_mesh_collection%get_global_mesh( global_mesh_names(i) )
    call global_mesh%get_target_mesh_names(target_mesh_names)
    if (allocated(target_mesh_names)) then
      do j=1, size(target_mesh_names)

        ! Only add a global mesh map if the required
        ! target mesh has been read in.
        if (global_mesh_collection%check_for(target_mesh_names(j))) then
          call file_handler%read_map( global_mesh_names(i), &
                                      target_mesh_names(j), &
                                      mesh_map )
          target_mesh => &
              global_mesh_collection%get_global_mesh( target_mesh_names(j) )
          call global_mesh%add_global_mesh_map( target_mesh, mesh_map )
        end if

      end do
    end if
  end do

  call file_handler%file_close()

  return
end subroutine assign_global_intergrid_maps


!>==============================================================================
!> @brief Sets common partition parameters to be applied to global meshes.
!>        (private subroutine)
!>
!> @param[in]   total_ranks  Total number of MPI ranks in this job
!> @param[out]  xproc              Number of ranks in mesh panel x-direction
!> @param[out]  yproc              Number of ranks in mesh panel y-direction
!> @param[out]  max_stencil_depth  Maximum depth of cells outside the base cell
!>                                 of stencil.
!> @param[out]  partitioner_ptr    Mesh partitioning strategy
!>==============================================================================
subroutine set_partition_parameters( total_ranks,       &
                                     xproc, yproc,      &
                                     max_stencil_depth, &
                                     partitioner_ptr )

  use base_mesh_config_mod,       only: geometry, &
                                        geometry_spherical
  use partitioning_config_mod,    only: panel_decomposition,        &
                                        panel_xproc, panel_yproc,   &
                                        PANEL_DECOMPOSITION_AUTO,   &
                                        PANEL_DECOMPOSITION_ROW,    &
                                        PANEL_DECOMPOSITION_COLUMN, &
                                        PANEL_DECOMPOSITION_CUSTOM
  use mixing_config_mod,          only: smagorinsky
  use partition_mod,              only: partitioner_cubedsphere_serial, &
                                        partitioner_cubedsphere,        &
                                        partitioner_planar
  use subgrid_config_mod,         only: dep_pt_stencil_extent, &
                                        rho_approximation_stencil_extent
  use transport_config_mod,       only: scheme,               &
                                        scheme_horz_cosmic,   &
                                        scheme_yz_bip_cosmic, &
                                        scheme_cosmic_3D,     &
                                        operators,            &
                                        operators_fv,         &
                                        fv_flux_order,        &
                                        fv_advective_order

  implicit none

  integer(i_def), intent(in)  :: total_ranks
  integer(i_def), intent(out) :: xproc
  integer(i_def), intent(out) :: yproc
  integer(i_def), intent(out) :: max_stencil_depth

  procedure(partitioner_interface), &
                  intent(out), pointer :: partitioner_ptr

  ! Locals
  integer(i_def) :: ranks_per_panel
  integer(i_def) :: start_factor
  integer(i_def) :: end_factor
  integer(i_def) :: fact_count
  integer(i_def) :: max_fv_stencil
  logical(l_def) :: found_factors

  character(str_def) :: domain_desc

  integer(i_def), parameter :: max_factor_iters = 10000

  partitioner_ptr => null()

  ! 1.0 Setup the partitioning strategy
  !===================================================================
  if (geometry == geometry_spherical) then

    if (total_ranks == 1 .or. mod(total_ranks,6) == 0) then

      ranks_per_panel = total_ranks/6
      domain_desc = "6x"

      if (total_ranks == 1) then
        ! Serial run job

        ranks_per_panel = 1
        partitioner_ptr => partitioner_cubedsphere_serial
        call log_event( "Using serial cubed sphere partitioner", &
                        LOG_LEVEL_INFO )

        if (scheme == scheme_horz_cosmic) then
          call log_event( "For Cosmic the total number of processors "//&
                          "must be >1 for a cubed-sphere domain.",      &
                          LOG_LEVEL_ERROR )
        end if

      else
        ! Paralled run job
        partitioner_ptr => partitioner_cubedsphere
        call log_event( "Using parallel cubed sphere partitioner", &
                        LOG_LEVEL_INFO )
      end if

    else
      call log_event( "Total number of processors must be 1 (serial) "//&
                      "or a multiple of 6 for a cubed-sphere domain.",  &
                      LOG_LEVEL_ERROR )
    end if

  else ! Planar/LAM mesh

    ranks_per_panel = total_ranks
    domain_desc = ""

    partitioner_ptr => partitioner_planar
    call log_event( "Using planar mesh partitioner ", &
                    LOG_LEVEL_INFO )
  end if


  ! 2.0 Setup Panel decomposition
  !===================================================================
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

  if (total_ranks > 1) then
    write(log_scratch_space, '(2(A,I0))' ) &
        'Panel decomposition: ', xproc,'x', yproc
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  end if

  ! 3.0 Determine max_stencil_depth
  !===================================================================
  max_stencil_depth = 1
  ! Smagorinsky (or boundary layers) appears to need larger haloes
  if (smagorinsky) max_stencil_depth = max(max_stencil_depth,2)
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

  return
end subroutine set_partition_parameters


!========================================================================
!> @brief Generates the 3D-mesh partitions required by the model
!>        configuration. (private subroutine)
!>
!> @details The extrusion types are setup for the required configuration
!>          before 3D-mesh partitions are instantiated.
!>
!> @param[in] local_rank   The local process number
!> @param[in] total_ranks  Total number of processes to run model
!> @param[in] xproc              Number of ranks in mesh panel x-direction
!> @param[in] yproc              Number of ranks in mesh panel y-direction
!> @param[in] max_stencil_depth  Maximum depth of cells outside the base cell
!>                               of stencil.
!> @param[in] partitioner_ptr    Mesh partitioning strategy
!> @param[in] create_2d_mesh           Create 2D mesh based on prime mesh
!> @param[in] create_shifted_mesh      Create shifted mesh based on prime mesh
!> @param[in] create_double_level_mesh Create double-level mesh based on prime mesh
!> @param[in] create_multigrid_meshes  Create meshes to support multigrid
!========================================================================
subroutine create_3D_mesh_partitions( local_rank, total_ranks,  &
                                      xproc, yproc,             &
                                      max_stencil_depth,        &
                                      partitioner_ptr,          &
                                      create_2d_mesh,           &
                                      create_shifted_mesh,      &
                                      create_double_level_mesh, &
                                      create_multigrid_meshes )

  use base_mesh_config_mod,    only: prime_mesh_name
  use extrusion_config_mod,    only: domain_top
  use formulation_config_mod,  only: l_multigrid
  use gungho_extrusion_mod,    only: create_extrusion,         &
                                     create_shifted_extrusion, &
                                     create_double_level_extrusion
  use mesh_collection_mod,     only: mesh_collection
  use multigrid_config_mod,    only: chain_mesh_tags

  implicit none

  integer(i_def), intent(in) :: local_rank
  integer(i_def), intent(in) :: total_ranks

  integer(i_def), intent(in) :: xproc
  integer(i_def), intent(in) :: yproc
  integer(i_def), intent(in) :: max_stencil_depth

  procedure(partitioner_interface), intent(in), pointer :: partitioner_ptr

  logical(l_def), intent(in) :: create_2d_mesh
  logical(l_def), intent(in) :: create_shifted_mesh
  logical(l_def), intent(in) :: create_double_level_mesh
  logical(l_def), intent(in) :: create_multigrid_meshes

  class(extrusion_type), allocatable :: extrusion
  class(extrusion_type), allocatable :: extrusion_shifted
  class(extrusion_type), allocatable :: extrusion_double
  type(uniform_extrusion_type)       :: extrusion_2d

  character(str_def) :: mesh_name
  character(str_def) :: source_mesh_name
  character(str_def) :: target_mesh_name

  integer(i_native) :: i

  integer(i_def), parameter :: one_layer = 1_i_def
  real(r_def),    parameter :: atmos_bottom = 0.0_r_def

  type(mesh_type), pointer :: source_mesh => null()
  type(mesh_type), pointer :: target_mesh => null()


  ! 1.0 Prime Mesh
  !===================================================================
  ! 1.1 Always generate the prime 3D mesh
  allocate( extrusion, &
            source=create_extrusion() )

  call create_3d_mesh_partition( local_rank, total_ranks, &
                                 xproc, yproc,            &
                                 max_stencil_depth,       &
                                 partitioner_ptr,         &
                                 prime_mesh_name,         &
                                 extrusion )

  ! 2.0 Generate addition 3d mesh partitions based on the prime mesh
  ! NOTE: This includes 2D meshes as they are currently impplemented
  !       as a 3D mesh of 1-layer thick.
  if (create_2d_mesh) then
    extrusion_2d = uniform_extrusion_type( atmos_bottom, &
                                           domain_top,   &
                                           one_layer )

    mesh_name = trim(prime_mesh_name)//'_2d'

    call create_3d_mesh_partition( local_rank, total_ranks, &
                                   xproc, yproc,            &
                                   max_stencil_depth,       &
                                   partitioner_ptr,         &
                                   prime_mesh_name,         &
                                   extrusion_2d,            &
                                   mesh_name=mesh_name )
  end if

  ! 3.0 Generate addition shifted 3d mesh partition
  !     based on the prime mesh.
  if (create_shifted_mesh) then

    allocate( extrusion_shifted, &
            source=create_shifted_extrusion(extrusion) )

    mesh_name = trim(prime_mesh_name)//'_shifted'

    call create_3d_mesh_partition( local_rank, total_ranks, &
                                   xproc, yproc,            &
                                   max_stencil_depth,       &
                                   partitioner_ptr,         &
                                   prime_mesh_name,         &
                                   extrusion_shifted,       &
                                   mesh_name=mesh_name )
  end if

  ! 4.0 Generate addition double-level 3d mesh partition
  !     based on the prime mesh.
  if (create_double_level_mesh) then

    allocate( extrusion_double, &
              source=create_double_level_extrusion(extrusion) )

    mesh_name = trim(prime_mesh_name)//'_double'

    call create_3d_mesh_partition( local_rank, total_ranks, &
                                   xproc, yproc,            &
                                   max_stencil_depth,       &
                                   partitioner_ptr,         &
                                   prime_mesh_name,         &
                                   extrusion_double,        &
                                   mesh_name=mesh_name )
  end if

  ! 5.0 Generate meshes required by any other schemes
  !     in the model run configuration.
  !===================================================================
  ! 5.0 Dynamics Multigrid
  if (create_multigrid_meshes) then

    do i=1, size(chain_mesh_tags)
      call create_3d_mesh_partition( local_rank, total_ranks, &
                                     xproc, yproc,            &
                                     max_stencil_depth,       &
                                     partitioner_ptr,         &
                                     chain_mesh_tags(i),      &
                                     extrusion )

      mesh_name = trim(chain_mesh_tags(i))//'_2d'
      call create_3d_mesh_partition( local_rank, total_ranks, &
                                     xproc, yproc,            &
                                     max_stencil_depth,       &
                                     partitioner_ptr,         &
                                     chain_mesh_tags(i),      &
                                     extrusion_2d,            &
                                     mesh_name=mesh_name )
    end do

    do i=1, size(chain_mesh_tags)-1
      source_mesh => mesh_collection%get_mesh(chain_mesh_tags(i))
      target_mesh => mesh_collection%get_mesh(chain_mesh_tags(i+1))

      call source_mesh%add_mesh_map(target_mesh)
      call target_mesh%add_mesh_map(source_mesh)

      source_mesh_name = trim(chain_mesh_tags(i))//'_2d'
      target_mesh_name = trim(chain_mesh_tags(i+1))//'_2d'
      source_mesh => mesh_collection%get_mesh(source_mesh_name)
      target_mesh => mesh_collection%get_mesh(target_mesh_name)

      call source_mesh%add_mesh_map(target_mesh)
      call target_mesh%add_mesh_map(source_mesh)

    end do

  end if

  if (allocated(extrusion))         deallocate( extrusion )
  if (allocated(extrusion_shifted)) deallocate( extrusion_shifted )
  if (allocated(extrusion_double))  deallocate( extrusion_double )

  return
end subroutine create_3D_mesh_partitions


!========================================================================
!> @brief   Generates a single 3D-mesh partition (private subroutine).
!> @details Instantiates a 3d-mesh partition and adds it to the model's
!>          mesh collection. Multiple meshes may be generated in the model
!>          based on the same global mesh but with differing extrusions.
!>
!> @param[in] local_rank            The local process number
!> @param[in] total_ranks           Total number of processes to run model
!> @param[in] xproc                 Number of ranks in mesh panel x-direction
!> @param[in] yproc                 Number of ranks in mesh panel y-direction
!> @param[in] max_stencil_depth     Maximum depth of cells outside the base cell
!>                                  of stencil.
!> @param[in] partitioner_ptr       Mesh partitioning strategy
!> @param[in] extrusion             Extrusion type for this 3D-mesh
!> @param[in] base_global_mesh_name Name of global base mesh
!> @param[in] mesh_name             Optional name of local 3D-mesh,
!>                                  defaults to name of base global mesh
!========================================================================
subroutine create_3D_mesh_partition( local_rank, total_ranks, &
                                     xproc, yproc,            &
                                     max_stencil_depth,       &
                                     partitioner_ptr,         &
                                     base_global_mesh_name,   &
                                     extrusion, mesh_name )

  use global_mesh_collection_mod, only: global_mesh_collection
  use partition_mod,              only: partitioner_cubedsphere_serial, &
                                        partitioner_cubedsphere,        &
                                        partitioner_planar
  use mesh_collection_mod,        only: mesh_collection

  implicit none

  integer(i_def), intent(in) :: local_rank
  integer(i_def), intent(in) :: total_ranks
  integer(i_def), intent(in) :: xproc
  integer(i_def), intent(in) :: yproc
  integer(i_def), intent(in) :: max_stencil_depth

  procedure(partitioner_interface), intent(in), pointer :: partitioner_ptr

  character(str_def),    intent(in) :: base_global_mesh_name
  class(extrusion_type), intent(in) :: extrusion
  character(str_def),    intent(in), &
                         optional   :: mesh_name

  type(global_mesh_type), pointer :: base_global_mesh => null()
  type(partition_type) :: partition

  type(mesh_type)    :: mesh
  integer(i_def)     :: mesh_id
  character(str_def) :: name

  if (.not. present(mesh_name)) then
    name = base_global_mesh_name
  else
    name = mesh_name
  end if

  ! 1.0 Check 3D-mesh hasn't already been created
  !===============================================
  if ( mesh_collection%check_for(name) ) return

  ! 2.0 Create the 3D-mesh
  !===============================================
  base_global_mesh => global_mesh_collection % &
                      get_global_mesh( base_global_mesh_name )

  partition = partition_type( base_global_mesh,  &
                              partitioner_ptr,   &
                              xproc, yproc,      &
                              max_stencil_depth, &
                              local_rank, total_ranks )

  mesh = mesh_type( base_global_mesh, &
                    partition,        &
                    extrusion,        &
                    mesh_name=name )

  mesh_id = mesh_collection%add_new_mesh( mesh )
  call mesh%clear()

  ! 3.0 Report on mesh creation
  !===============================================
  write(log_scratch_space,'(A,I0,A)')                 &
      '   ... "'//trim(name)//'"(id:', mesh_id,') '// &
      'based on global mesh "'//trim(base_global_mesh_name)//'"'

  if (mesh_id /= imdi) then
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  else
    write(log_scratch_space,'(A,I0,A)') &
        trim(log_scratch_space)//' (FAILED)'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  return
end subroutine create_3D_mesh_partition


!========================================================================
!> @brief   Assigns intergrid maps to 3D-mesh partitions
!>          (private subroutine).
!> @details Adds local ID integrid mappings 3D-mesh parititons.
!>          Local integrid maps are assign to to source (source-> target)
!>          and target (target -> source) meshes. A number of
!>          of assumptions are made.
!>
!>          *  The base global meshes of source/target 3D-mesh
!>             are present in the global_mesh_collection.
!>          *  Integrid maps for the base global meshes where
!>             where read in and assigned.
!>          *  Partitioned 3D-mesh are present in the
!>             mesh_collection.
!>
!> @param[in] source_mesh_name   Name of source 3D mesh partition
!> @param[in] target_mesh_name   Name of target 3D mesh partition
!========================================================================
subroutine add_intergrid_maps( source_mesh_name, &
                               target_mesh_name )

  use mesh_collection_mod, only: mesh_collection

  implicit none

  character(str_def), intent(in) :: source_mesh_name
  character(str_def), intent(in) :: target_mesh_name

  type(mesh_type), pointer :: source_mesh => null()
  type(mesh_type), pointer :: target_mesh => null()


  ! Now add in any mesh maps required by multigrid
  source_mesh => mesh_collection % get_mesh( source_mesh_name )
  target_mesh => mesh_collection % get_mesh( target_mesh_name )

  if ( associated(source_mesh) .and. &
       associated(target_mesh) ) then

    ! Mesh tag names may be different but "could point to the same mesh
    ! So check the ids are not the same
    if (source_mesh%get_id() == target_mesh%get_id()) then
      write(log_scratch_space,'(A)')                  &
          'Unable to create intergrid map: Source('// &
          trim(source_mesh_name)//' and target('//    &
          trim(target_mesh_name)//') mesh ids are the same'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR)
    end if

    call source_mesh % add_mesh_map (target_mesh)
    call target_mesh % add_mesh_map (source_mesh)
    write(log_scratch_space,'(A,I0,A)')     &
        'Adding intergrid map "'//          &
         trim(source_mesh_name)//'"<-->"'// &
         trim(target_mesh_name)//'"'
    call log_event( log_scratch_space, LOG_LEVEL_INFO)
  else
    write(log_scratch_space,'(A,I0,A)')          &
        'Unable to create mesh map between "'//  &
        trim(source_mesh_name)//'"-"'//          &
        trim(target_mesh_name)//'"'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR)
  end if

  nullify(source_mesh)
  nullify(target_mesh)

  return
end subroutine add_intergrid_maps


!========================================================================
!> @brief Loads list of global meshes from configuration mesh input
!>        file (private subroutine).
!>
!> @param[in] mesh_names[:]   Array of requested mesh names to load
!>                            from the configuration mesh input
!>                            file.
!> @param[in] n_panels        Number of panel domains in global mesh
!========================================================================
subroutine load_mesh(mesh_names, n_panels)

  use base_mesh_config_mod,       only: filename
  use global_mesh_collection_mod, only: global_mesh_collection

  implicit none

  character(str_def), intent(in) :: mesh_names(:)
  integer(i_def),     intent(in) :: n_panels

  type(ugrid_mesh_data_type) :: ugrid_mesh_data
  type (global_mesh_type)    :: global_mesh

  integer(i_def) :: i

  do i=1, size(mesh_names)
    if (.not. global_mesh_collection%check_for(mesh_names(i))) then

      call ugrid_mesh_data%read_from_file(trim(filename), mesh_names(i))
      global_mesh = global_mesh_type( ugrid_mesh_data, n_panels )
      call ugrid_mesh_data%clear()

      call global_mesh_collection%add_new_global_mesh ( global_mesh )

    end if
  end do

  return
end subroutine load_mesh


!========================================================================
!> @brief Finalises the mesh_collection
!========================================================================
subroutine final_mesh()

  use mesh_collection_mod, only: mesh_collection

  implicit none

  if (allocated(mesh_collection)) then
    call mesh_collection%clear()
    deallocate(mesh_collection)
  end if

  return
end subroutine final_mesh

end module create_mesh_mod
