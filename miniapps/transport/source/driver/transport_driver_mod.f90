!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>
!> @brief Drives the execution of the transport miniapp.
!>
module transport_driver_mod

  use base_mesh_config_mod,           only: prime_mesh_name
  use checksum_alg_mod,               only: checksum_alg
  use clock_mod,                      only: clock_type
  use configuration_mod,              only: final_configuration
  use constants_mod,                  only: i_def, r_def, l_def, i_native, &
                                            tiny_eps, l_def
  use convert_to_upper_mod,           only: convert_to_upper
  use create_fem_mod,                 only: init_fem
  use create_mesh_mod,                only: init_mesh
  use derived_config_mod,             only: set_derived_config
  use field_mod,                      only: field_type
  use global_mesh_collection_mod,     only: global_mesh_collection, &
                                            global_mesh_collection_type
  use mesh_collection_mod,            only: mesh_collection
  use init_clock_mod,                 only: initialise_clock
  use init_transport_mod,             only: init_transport
  use transport_mod,                  only: transport_load_configuration, &
                                            program_name
  use lfric_xios_io_mod,              only: initialise_xios
  use diagnostics_io_mod,             only: write_scalar_diagnostic, &
                                            write_vector_diagnostic
  use diagnostics_calc_mod,           only: write_density_diagnostic
  use local_mesh_collection_mod,      only: local_mesh_collection, &
                                            local_mesh_collection_type
  use log_mod,                        only: log_event,                        &
                                            log_set_level,                    &
                                            log_scratch_space,                &
                                            initialise_logging,               &
                                            finalise_logging,                 &
                                            LOG_LEVEL_ALWAYS,                 &
                                            LOG_LEVEL_ERROR,                  &
                                            LOG_LEVEL_WARNING,                &
                                            LOG_LEVEL_INFO,                   &
                                            LOG_LEVEL_DEBUG,                  &
                                            LOG_LEVEL_TRACE
  use mesh_collection_mod,            only: mesh_collection
  use io_config_mod,                  only: diagnostic_frequency,             &
                                            nodal_output_on_w3,               &
                                            write_diag,                       &
                                            use_xios_io,                      &
                                            subroutine_timers
  use timer_mod,                      only: init_timer, timer, output_timer
  use mpi_mod,                        only: store_comm,    &
                                            get_comm_size, &
                                            get_comm_rank
  use transport_config_mod,           only: key_from_scheme,            &
                                            scheme,                     &
                                            scheme_yz_bip_cosmic,       &
                                            scheme_cosmic_3D,           &
                                            scheme_horz_cosmic,         &
                                            scheme_method_of_lines,     &
                                            rho_splitting,              &
                                            rho_splitting_none,         &
                                            theta_splitting,            &
                                            theta_splitting_none,       &
                                            cfl_mol_3d_stab
  use mass_conservation_alg_mod,      only: mass_conservation
  use yz_bip_cosmic_alg_mod,          only: yz_bip_cosmic_step
  use cusph_cosmic_transport_alg_mod, only: set_winds,      &
                                            cusph_cosmic_transport_step
  use cosmic_threed_alg_mod,          only: cosmic_threed_transport_step
  use calc_dep_pts_alg_mod,           only: calc_dep_pts
  use density_inc_update_alg_mod,     only: density_inc_update_alg
  use runtime_constants_mod,          only: get_detj_at_w2,                   &
                                            get_detj_at_w2_shifted,           &
                                            get_cell_orientation,             &
                                            get_cell_orientation_shifted
  use yaxt,                           only: xt_initialize, xt_finalize
  use xios,                           only: xios_context_finalize, &
                                            xios_update_calendar
  use rk_transport_rho_mod,           only: rk_transport_rho_final
  use rk_transport_theta_mod,         only: rk_transport_theta_final
  use runge_kutta_init_mod,           only: runge_kutta_init, &
                                            runge_kutta_final
  use split_mol_cosmic_alg_mod,       only: split_transport_rho_step
  use diagnostic_alg_mod,             only: l2norm_diff_2fields
  use advection_alg_mod,              only: advection_alg_init
  use transport_general_mod,          only: transport_general
  use initial_wind_config_mod,        only: profile,                      &
                                            profile_curl_free_reversible, &
                                            profile_div_free_reversible,  &
                                            profile_xy_NL_case_1,         &
                                            profile_yz_NL_case_1,         &
                                            profile_NL_case_1,            &
                                            profile_NL_case_2,            &
                                            profile_NL_case_3,            &
                                            profile_NL_case_4,            &
                                            profile_hadley_like_dcmip


  implicit none

  private

  public :: initialise_transport, run_transport, finalise_transport

  ! Prognostic fields
  type(field_type) :: wind_n
  type(field_type) :: density
  type(field_type) :: theta
  type(field_type) :: dep_pts_x
  type(field_type) :: dep_pts_y
  type(field_type) :: dep_pts_z

  ! Shifted Prognostic fields
  type(field_type) :: wind_shifted
  type(field_type) :: density_shifted
  type(field_type) :: increment
  type(field_type) :: wind_divergence
  type(field_type), pointer :: detj_at_w2 => null()
  type(field_type), pointer :: detj_at_w2_shifted => null()
  type(field_type), pointer :: cell_orientation => null()
  type(field_type), pointer :: cell_orientation_shifted => null()

  ! Coordinate field
  type(field_type), target, dimension(3) :: chi_xyz
  type(field_type), target, dimension(3) :: chi_sph
  type(field_type), target               :: panel_id
  type(field_type), target, dimension(3) :: shifted_chi_xyz
  type(field_type), target, dimension(3) :: shifted_chi_sph

  class(clock_type), allocatable :: clock

  integer(i_def) :: mesh_id
  integer(i_def) :: twod_mesh_id
  integer(i_def) :: shifted_mesh_id

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Sets up required state in preparation for run.
  !! @param[in] filename Configuration namelist file
  subroutine initialise_transport( filename, model_communicator )

    use logging_config_mod, only: run_log_level,          &
                                  key_from_run_log_level, &
                                  RUN_LOG_LEVEL_ERROR,    &
                                  RUN_LOG_LEVEL_INFO,     &
                                  RUN_LOG_LEVEL_DEBUG,    &
                                  RUN_LOG_LEVEL_TRACE,    &
                                  RUN_LOG_LEVEL_WARNING

    implicit none

    character(*),      intent(in) :: filename
    integer(i_native), intent(in) :: model_communicator

    character(len=*), parameter :: xios_ctx  = "transport"

    integer(i_def)    :: total_ranks, local_rank
    integer(i_native) :: log_level

    integer(i_def), allocatable :: multigrid_mesh_ids(:)
    integer(i_def), allocatable :: multigrid_2d_mesh_ids(:)

    ! Store the MPI communicator for later use
    call store_comm( model_communicator )

    call transport_load_configuration( filename )

    call initialise_clock( clock )

    ! Initialise YAXT
    call xt_initialize( model_communicator )

    ! Get the rank information
    total_ranks = get_comm_size()
    local_rank  = get_comm_rank()

    call initialise_logging( local_rank, total_ranks, program_name )

    select case (run_log_level)
    case( RUN_LOG_LEVEL_ERROR )
      log_level = LOG_LEVEL_ERROR
    case( RUN_LOG_LEVEL_WARNING )
      log_level = LOG_LEVEL_WARNING
    case( RUN_LOG_LEVEL_INFO )
      log_level = LOG_LEVEL_INFO
    case( RUN_LOG_LEVEL_DEBUG )
      log_level = LOG_LEVEL_DEBUG
    case( RUN_LOG_LEVEL_TRACE )
      log_level = LOG_LEVEL_TRACE
    end select

    call log_set_level( log_level )

    write(log_scratch_space,'(A)')                             &
       'Runtime message logging severity set to log level: '// &
       convert_to_upper(key_from_run_log_level(run_log_level))
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    call log_event( program_name//': Runtime default precision set as:', LOG_LEVEL_ALWAYS )
    write(log_scratch_space, '(I1)') kind(1.0_r_def)
    call log_event( program_name//':     r_def kind = '//log_scratch_space , LOG_LEVEL_ALWAYS )
    write(log_scratch_space, '(I1)') kind(1_i_def)
    call log_event( program_name//':     i_def kind = '//log_scratch_space , LOG_LEVEL_ALWAYS )

    call set_derived_config( .true. )

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------
    call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )
    if ( subroutine_timers ) then
      call init_timer()
      call timer( program_name )
    end if

    ! Mesh initialisation
    allocate( global_mesh_collection, &
              source = global_mesh_collection_type() )

    allocate( local_mesh_collection, &
              source = local_mesh_collection_type() )

    ! Create the mesh
    call init_mesh( local_rank, total_ranks, mesh_id,      &
                    twod_mesh_id=twod_mesh_id,             &
                    shifted_mesh_id=shifted_mesh_id,       &
                    multigrid_mesh_ids=multigrid_mesh_ids, &
                    multigrid_2d_mesh_ids=multigrid_2d_mesh_ids )

    ! FEM initialisation
    call init_fem( mesh_id, chi_xyz, chi_sph, panel_id,    &
                   shifted_mesh_id=shifted_mesh_id,        &
                   shifted_chi_xyz=shifted_chi_xyz,        &
                   shifted_chi_sph=shifted_chi_sph,        &
                   multigrid_mesh_ids=multigrid_mesh_ids,  &
                   multigrid_2d_mesh_ids=multigrid_2d_mesh_ids )

    ! Transport initialisation
    call init_transport( mesh_id, twod_mesh_id, chi_xyz, chi_sph, panel_id, &
                         shifted_mesh_id, shifted_chi_xyz, shifted_chi_sph, &
                         wind_n, density, theta, dep_pts_x, dep_pts_y,      &
                         dep_pts_z, increment, wind_divergence,             &
                         wind_shifted, density_shifted )

    ! Full global meshes no longer required, so reclaim
    ! the memory from global_mesh_collection
    if (allocated(global_mesh_collection)) deallocate(global_mesh_collection)

    ! Initialise dependencies
    ! scheme_method_of_lines                  -> runge_kutta_init
    ! rho_splitting /= rho_splitting_none   or
    ! theta_splitting /= theta_splitting_none -> advection_alg_init
    !
    if ( scheme == scheme_method_of_lines    .or.  &
         rho_splitting /= rho_splitting_none .or.  &
         theta_splitting /= theta_splitting_none   ) then
      call runge_kutta_init()
      call advection_alg_init( mesh_id )
    end if

    ! Calculate det(J) at W2 dofs for chi and shifted_chi fields.
    ! The calculation of det(J) for the shifted_chi field is done in preparation
    ! for Ticket #1608.
    detj_at_w2 => get_detj_at_w2()
    detj_at_w2_shifted => get_detj_at_w2_shifted()

    ! Set cell cell_orientation (and shifted) if using splitting scheme.
    ! Cell orientation is used by Cosmic and splitting schemes
    !
    if (scheme /= scheme_method_of_lines    .or. &
        rho_splitting /= rho_splitting_none .or. &
        theta_splitting /= theta_splitting_none  ) then
      cell_orientation => get_cell_orientation()
      cell_orientation_shifted => get_cell_orientation_shifted()
    end if

    if ( use_xios_io ) then
      call initialise_xios( xios_ctx,     &
                            model_communicator, &
                            clock,        &
                            mesh_id,      &
                            twod_mesh_id, &
                            chi_xyz )
    end if

    ! Output initial conditions
    if (clock%is_initialisation()) then

      if ( write_diag ) then

        if ( use_xios_io ) then

          ! Need to ensure calendar is initialised here as XIOS has no concept of timestep 0
          call xios_update_calendar(1)

        end if

        call write_vector_diagnostic( 'wind', wind_n, clock, &
                                      mesh_id, nodal_output_on_w3 )
        call write_scalar_diagnostic( 'density', density, clock, &
                                      mesh_id, nodal_output_on_w3 )
        call write_scalar_diagnostic( 'theta', theta, clock, &
                                      mesh_id, nodal_output_on_w3 )

     end if

    end if

  end subroutine initialise_transport

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Performs time stepping.
  !>
  subroutine run_transport()

    implicit none

    type(field_type)    :: density_t0, density_inc, density_n
    type(field_type)    :: theta_t0, theta_inc, theta_n
    real(r_def)         :: err_rho, err_theta
    logical(l_def)      :: conservative_form, shifted_grid
    logical(l_def)      :: time_varying_wind


    call density_t0%initialise( vector_space = density%get_function_space() )
    call density_inc%initialise( vector_space = density%get_function_space() )
    call density_n%initialise( vector_space = density%get_function_space() )

    call theta_t0%initialise( vector_space = theta%get_function_space() )
    call theta_inc%initialise( vector_space = theta%get_function_space() )
    call theta_n%initialise( vector_space = theta%get_function_space() )

    call density%copy_field(density_t0)
    call theta%copy_field(theta_t0)

    call log_event( 'Running '//program_name//' ...', LOG_LEVEL_ALWAYS )
    call log_event( program_name//': Miniapp will run with default precision set as:', LOG_LEVEL_INFO )
    write(log_scratch_space, '(I1)') kind(1.0_r_def)
    call log_event( program_name//':        r_def kind = '//log_scratch_space , LOG_LEVEL_INFO )
    write(log_scratch_space, '(I1)') kind(1_i_def)
    call log_event( program_name//':        i_def kind = '//log_scratch_space , LOG_LEVEL_INFO )

    if ( theta_splitting == theta_splitting_none ) then
      ! Currently there is no other theta transport option in the tranport miniapp except splitting
      write( log_scratch_space, '(A)' ) ' WARNING == Invalid option for theta advection for the miniapp transport'
      call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      write( log_scratch_space, '(A)' ) ' In this case theta is not advected and remains at the initial state '
      call log_event( log_scratch_space, LOG_LEVEL_WARNING )
    end if

    ! Determine if wind is time-varying or not
    if ((profile == profile_xy_NL_case_1)         .or. &
        (profile == profile_yz_NL_case_1)         .or. &
        (profile == profile_NL_case_1)            .or. &
        (profile == profile_NL_case_2)            .or. &
        (profile == profile_NL_case_3)            .or. &
        (profile == profile_NL_case_4)            .or. &
        (profile == profile_hadley_like_dcmip)    .or. &
        (profile == profile_div_free_reversible)  .or. &
        (profile == profile_curl_free_reversible)) then
        time_varying_wind = .true.
      else
        time_varying_wind = .false.
      end if

    ! Initialise winds before first time step
    call set_winds( wind_n, mesh_id, clock%get_step() )

    !--------------------------------------------------------------------------
    ! Model step
    !--------------------------------------------------------------------------
    do while (clock%tick())

      write(log_scratch_space, '("/", A, "\ ")') repeat('*', 76)
      call log_event( log_scratch_space, LOG_LEVEL_TRACE )
      write( log_scratch_space, '(A,I0)' ) &
        'Start of timestep ', clock%get_step()
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      ! Update XIOS calendar if we are using it for diagnostic output or
      ! checkpoint
      if ( use_xios_io ) then
        call log_event( "Transport: Updating XIOS timestep", LOG_LEVEL_INFO )
        call xios_update_calendar( clock%get_step() )
      end if

      ! Only update the wind for time varying prescribed profiles
      if (time_varying_wind) then
        call set_winds( wind_n, mesh_id, clock%get_step() )
      end if

      if (scheme /= scheme_method_of_lines    .or. &
          rho_splitting /= rho_splitting_none .or. &
          theta_splitting /= theta_splitting_none  ) then
        ! Calculate departure points.
        ! Here the wind is assumed to be the same at timestep n and timestep n+1
        call calc_dep_pts( dep_pts_x, dep_pts_y, dep_pts_z, wind_divergence,    &
                           wind_n, wind_n, detj_at_w2, cell_orientation )
      end if

      if ( subroutine_timers ) call timer( 'transport step' )

      call density%copy_field(density_n)

      ! Transport of rho field
      if ( rho_splitting == rho_splitting_none ) then
        ! No splitting, so use a '3D' implementation
        select case( scheme )
          case ( scheme_yz_bip_cosmic )
            call yz_bip_cosmic_step( increment, density, dep_pts_y, dep_pts_z, &
                                                                   detj_at_w2  )
            call density_inc_update_alg(density, increment)

          case ( scheme_horz_cosmic )
            call cusph_cosmic_transport_step( increment, density, dep_pts_x,    &
                                        dep_pts_y, detj_at_w2, cell_orientation )
            call density_inc_update_alg(density, increment)

          case ( scheme_cosmic_3D )
            call cosmic_threed_transport_step( increment, density, dep_pts_x,   &
                             dep_pts_y, dep_pts_z, detj_at_w2, cell_orientation )
            call density_inc_update_alg(density, increment)

          case ( scheme_method_of_lines )

            ! shifted_grid = T/F (the field is on shifted/not grid)
            ! conservative_form=T/F (solving the conservative/advective form of transport)
            shifted_grid = .false.
            conservative_form = .true.

            call transport_general(wind_n, density, scheme_method_of_lines, &
                                   shifted_grid, conservative_form          )
            conservative_form = .false.
            call transport_general(wind_n, theta, scheme_method_of_lines,   &
                                  shifted_grid, conservative_form           )

          case default
            write(log_scratch_space, '(A, A)') &
               "Unrecognised transport scheme: ", key_from_scheme(scheme)
            call log_event( log_scratch_space, LOG_LEVEL_ERROR )
            stop
        end select
      else
        ! Use the split scheme for density and theta fields
        !-----------------------------------------------------
        ! Both theta and density are not on the shifted grid (shifted_grid = .false.)
        ! Density is in the conservative for (mconservative_form = .true.
        ! Theta is in advection form (conservative_form = .false.)

         shifted_grid = .false.
         conservative_form = .true.
         call transport_general(wind_n, density, rho_splitting, &
                                shifted_grid, conservative_form )
         conservative_form = .false.
         call transport_general(wind_n, theta, theta_splitting, &
                                shifted_grid, conservative_form )
      end if

      if ( subroutine_timers ) call timer( 'transport step' )

      call write_density_diagnostic( density, clock )
      call mass_conservation( clock%get_step(), density_n, density )

      call density%log_minmax( LOG_LEVEL_INFO, 'density' )
      call   theta%log_minmax( LOG_LEVEL_INFO, 'theta'   )
      call l2norm_diff_2fields(err_rho, density, density_t0 )
      call l2norm_diff_2fields(err_theta, theta, theta_t0 )

      write( log_scratch_space, '(A,E16.8)' ) ' L2(rho-rho0)/L2(rho0)       = ',err_rho
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
      write( log_scratch_space, '(A,E16.8)' ) ' L2(theta-theta0)/L2(theta0) = ',err_theta
      call log_event( log_scratch_space, LOG_LEVEL_INFO )


      write( log_scratch_space, '(A,I0)' ) 'End of timestep ', clock%get_step()
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
      write(log_scratch_space, '("\", A, "/ ")') repeat('*', 76)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      ! Output wind and density values.
      if ( (mod( clock%get_step(), diagnostic_frequency ) == 0) &
           .and. write_diag ) then

        call write_vector_diagnostic( 'wind', wind_n,                     &
                                      clock, mesh_id, nodal_output_on_w3)
        call write_scalar_diagnostic( 'density', density,                 &
                                      clock, mesh_id, nodal_output_on_w3)
        call write_scalar_diagnostic( 'wind_divergence', wind_divergence, &
                                      clock, mesh_id, nodal_output_on_w3)
        call write_scalar_diagnostic( 'theta', theta,                     &
                                      clock, mesh_id, nodal_output_on_w3)

      end if

    end do ! while clock%is_running()

    call l2norm_diff_2fields(err_rho, density, density_t0 )
    call l2norm_diff_2fields(err_theta, theta, theta_t0 )

    write( log_scratch_space, '(A,E16.8)' ) ' l2-norm error (rho)   @ END = ',err_rho
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(A,E16.8)' ) ' l2-norm error (theta) @ END = ',err_theta
    call log_event( log_scratch_space, LOG_LEVEL_INFO )


    if ( scheme == scheme_method_of_lines    .or.  &
         rho_splitting /= rho_splitting_none .or.  &
         theta_splitting /= theta_splitting_none   ) then
      call runge_kutta_final()
      call rk_transport_rho_final()
      call rk_transport_theta_final()
    end if

    nullify( detj_at_w2, detj_at_w2_shifted, cell_orientation, &
             cell_orientation_shifted )

  end subroutine run_transport

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Tidies up after a run.
  !>
  subroutine finalise_transport()

    implicit none

    !----------------------------------------------------------------------------
    ! Model finalise
    !----------------------------------------------------------------------------
    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Write checksums to file
    call checksum_alg( program_name, density, 'density',  &
                                     wind_n, 'wind',      &
                                     theta, 'theta'       )

    if ( subroutine_timers ) then
      call timer( program_name )
      call output_timer()
    end if

    if ( use_xios_io ) then
      call xios_context_finalize()
    end if

    ! Finalise namelist configurations
    call final_configuration()

    !----------------------------------------------------------------------------
    ! Driver layer finalise
    !----------------------------------------------------------------------------
    ! Finalise YAXT
    call xt_finalize()

    call log_event( program_name//' completed.', LOG_LEVEL_ALWAYS )

    ! Finalise the logging system
    call finalise_logging()

  end subroutine finalise_transport

end module transport_driver_mod
