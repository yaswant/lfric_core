!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Container for Gung Ho model run working data set inclduing methods
!> to initialise, copy and finalise the data set
!>
!> This module provides a type to hold all the model fields and methods to
!> initialise (create and read), copy and finalise (write and destroy) the
!> data contained within the type.
!>
module gungho_model_data_mod

  use mr_indices_mod,                     only : nummr
  use moist_dyn_mod,                      only : num_moist_factors
  use field_mod,                          only : field_type
  use field_parent_mod,                   only : write_interface
  use field_collection_mod,               only : field_collection_type
  use constants_mod,                      only : i_def, l_def, r_def
  use log_mod,                            only : log_event,       &
                                                 LOG_LEVEL_INFO,  &
                                                 LOG_LEVEL_ERROR, &
                                                 LOG_LEVEL_DEBUG, &
                                                 log_scratch_space
  use base_mesh_config_mod,               only : geometry,          &
                                                 geometry_spherical
  use mesh_mod,                           only : mesh_type
  use files_config_mod,                   only : checkpoint_stem_name
  use formulation_config_mod,             only : use_physics
  use initialization_config_mod,          only : init_option,                 &
                                                 init_option_analytic,        &
                                                 init_option_fd_start_dump,   &
                                                 init_option_checkpoint_dump, &
                                                 init_option_fe_start_dump,   &
                                                 ancil_option,                &
                                                 ancil_option_none,           &
                                                 ancil_option_start_dump,     &
                                                 ancil_option_fixed,          &
                                                 ancil_option_updating,       &
                                                 lbc_option,                  &
                                                 lbc_option_none,             &
                                                 lbc_option_analytic,         &
                                                 lbc_option_gungho_file,      &
                                                 lbc_option_um2lfric_file
  use io_config_mod,                      only : checkpoint_read,  &
                                                 checkpoint_write, &
                                                 write_dump
  use lfric_xios_read_mod,                only : read_checkpoint,  &
                                                 read_state
  use lfric_xios_write_mod,               only : write_checkpoint, &
                                                 write_state
  use boundaries_config_mod,              only : limited_area
  use create_lbcs_mod,                    only : create_lbc_fields
  use create_gungho_prognostics_mod,      only : create_gungho_prognostics
  use create_physics_prognostics_mod,     only : create_physics_prognostics
  use section_choice_config_mod,          only : cloud, cloud_none
  use map_fd_to_prognostics_alg_mod,      only : map_fd_to_prognostics
  use model_clock_mod,                    only : model_clock_type
  use gungho_init_prognostics_driver_mod, only : init_gungho_prognostics
  use init_gungho_lbcs_alg_mod,           only : init_lbcs_file_alg,    &
                                                 init_lbcs_analytic_alg
  use init_physics_prognostics_alg_mod,   only : init_physics_prognostics_alg
  use derived_config_mod,                 only : l_esm_couple
#ifdef COUPLED
  use coupler_mod,                      only : cpl_init_fields
#endif
  use moist_dyn_factors_alg_mod,        only : moist_dyn_factors_alg
  use init_fd_prognostics_mod,          only : init_fd_prognostics_dump
  use map_physics_fields_alg_mod,       only : map_physics_fields_alg
  use set_any_dof_alg_mod,              only : set_any_dof_alg
  use reference_element_mod,            only : T
#ifdef UM_PHYSICS
  use create_fd_prognostics_mod,          only : create_fd_prognostics
  use init_ancils_mod,                    only : create_fd_ancils
  use process_inputs_alg_mod,             only : process_inputs_alg
  use update_tstar_alg_mod,               only : update_tstar_alg
  use variable_fields_mod,                only : init_variable_fields
#endif
  use linked_list_mod,                    only : linked_list_type, &
                                                 linked_list_item_type
  use energy_correction_config_mod,       only : encorr_usage, &
                                                 encorr_usage_none
  use compute_total_energy_alg_mod,       only : compute_total_energy_alg
  use compute_total_mass_alg_mod,         only : compute_total_mass_alg
  use compute_column_integral_kernel_mod, &
                                          only : compute_column_integral_kernel_type
  use geometric_constants_mod,            only : get_height
  use planet_config_mod,                  only : radius
  use scalar_to_field_alg_mod,            only : scalar_to_field_alg
  use field_to_scalar_alg_mod,            only : field_to_scalar_alg

  implicit none

  private

  !> Holds the working data set for a model run and other working state.
  !>
  type :: model_data_type

    private

    !> Stores all the fields used by the model - the remaining collections store
    !> pointers to the data in the depository.
    type( field_collection_type ), public :: depository

    !> @name Fields needed to time-step the model.
    !> @{

    !> All the prognostic fields (except for field arrays: auxiliary prognostic)
    type( field_collection_type ), public   :: prognostic_fields
    !> All the diagnostic fields
    type( field_collection_type ), public   :: diagnostic_fields
    !> Fields that should be advected
    type( field_collection_type ), public   :: adv_fields_last_outer
    type( field_collection_type ), public   :: adv_fields_all_outer
    !> FD fields derived from FE fields for use in physics time-stepping schemes
    type( field_collection_type ), public   :: derived_fields
    !> LBC fields - lateral boundary conditions to run a limited area model
    type( field_collection_type ), public   :: lbc_fields
    !> Linked list of time axis objects used to update time-varying lbcs
    type( linked_list_type ),      public   :: lbc_times_list
    !> Fields owned by the radiation scheme
    type( field_collection_type ), public   :: radiation_fields
    !> Fields owned by the microphysics scheme
    type( field_collection_type ), public   :: microphysics_fields
    !> Fields owned by the electric scheme
    type( field_collection_type ), public   :: electric_fields
    !> Fields owned by the orographic drag schemes
    type( field_collection_type ), public   :: orography_fields
    !> Fields owned by the turbulence scheme
    type( field_collection_type ), public   :: turbulence_fields
    !> Fields owned by the convection schemes
    type( field_collection_type ), public   :: convection_fields
    !> Fields owned by the cloud schemes
    type( field_collection_type ), public   :: cloud_fields
    !> Fields owned by the surface exchange scheme
    type( field_collection_type ), public   :: surface_fields
    !> Fields owned by the soil hydrology scheme
    type( field_collection_type ), public   :: soil_fields
    !> Fields owned by the snow scheme
    type( field_collection_type ), public   :: snow_fields
    !> Fields owned by the chemistry schemes
    type( field_collection_type ), public   :: chemistry_fields
    !> Fields owned by the aerosol schemes
    type( field_collection_type ), public   :: aerosol_fields
    !> Fields owned by the stochastic physics schemes
    type( field_collection_type ), public   :: stph_fields
    !> Array of fields containing the moisture mixing ratios
    !>  (auxiliary prognostic)
    type( field_type ), allocatable, public :: mr(:)
    !> Array of fields containing the moist dynamics (auxiliary prognostic)
    type( field_type ), allocatable, public :: moist_dyn(:)
    !> Array of fields containing coupling data
    type( field_collection_type ), public :: cpl_snd
    type( field_collection_type ), public :: cpl_rcv
    !> @}

    !> FD fields used to read initial conditions from LFRic-Input files
    type( field_collection_type ), public   :: fd_fields

    !> Fields used to store data read in from ancillary files
    type( field_collection_type ), public   :: ancil_fields
    !> Linked list of time axis objects used to update time-varying ancils
    type( linked_list_type ),      public   :: ancil_times_list

    !> Fields for the tangent linear linearisation state
    type( field_collection_type ), public   :: ls_fields
    !> Array of linearisation fields containing the moisture mixing ratios
    type( field_type ), allocatable, public :: ls_mr(:)
    !> Array of linearisation fields containing the moist dynamics
    type( field_type ), allocatable, public :: ls_moist_dyn(:)
    !> Linked list of time axis objects used to update time-varying
    !> linearisation state
    type( linked_list_type ),      public   :: ls_times_list

    !> Rate of temperature adjustment for energy correction
    real( kind=r_def ), public :: temperature_correction_rate
    !> Total mass of dry atmosphere used for energy correction
    real( kind=r_def ), public :: total_dry_mass
    !> Total energy of moist atmosphere for calculating energy correction
    real( kind=r_def ), public :: total_energy
    !> Total energy of moist atmosphere at previous energy correction step
    real( kind=r_def ), public :: total_energy_previous

    contains

  end type model_data_type

  ! Set these to select how to initialize model prognostic fields
  integer(i_def) :: prognostic_init_choice, ancil_choice

  public model_data_type, create_model_data, finalise_model_data, &
         initialise_model_data, output_model_data

contains

  !> @brief Create the fields contained in model_data
  !> @param[inout] model_data The working data set for a model run
  !> @param[in]    mesh      The current 3d mesh
  !> @param[in]    twod_mesh The current 2d mesh
  !> @param[in]    aerosol_mesh      Aerosol 3d mesh
  !> @param[in]    aerosol_twod_mesh Aerosol 2d mesh
subroutine create_model_data( model_data,         &
                              mesh,               &
                              twod_mesh,          &
                              aerosol_mesh,       &
                              aerosol_twod_mesh,  &
                              model_clock         )


    implicit none

    type( model_data_type ), intent(inout)       :: model_data
    type( mesh_type ),       intent(in), pointer :: mesh
    type( mesh_type ),       intent(in), pointer :: twod_mesh
    type( mesh_type ),       intent(in), pointer :: aerosol_mesh
    type( mesh_type ),       intent(in), pointer :: aerosol_twod_mesh
    class(model_clock_type), intent(in)          :: model_clock

    !-------------------------------------------------------------------------
    ! Select how to initialize model prognostic fields
    !-------------------------------------------------------------------------

    ! This way of setting up the initialisation options is not ideal, but
    ! pragmatic for now and avoids extra namelist changes. It should be
    ! reviewed in the next round of driver layer refactoring

    ! Get the specified namelist options for prognostic initialisation
    prognostic_init_choice = init_option
    ancil_choice = ancil_option

    ! If checkpoint reading has been specified then override these options
    if (checkpoint_read) then
      prognostic_init_choice = init_option_checkpoint_dump
      if ( ancil_option /= ancil_option_updating ) then
        ancil_choice = ancil_option_none
      end if
    end if

    !-------------------------------------------------------------------------
    ! Instantiate the fields
    !-------------------------------------------------------------------------

    ! Field bundles - allocate the fields so thay can be cleared
    allocate(model_data%moist_dyn(num_moist_factors))
    allocate(model_data%mr(nummr))
    allocate(model_data%ls_moist_dyn(num_moist_factors))
    allocate(model_data%ls_mr(nummr))

    ! Create gungho prognostics and auxilliary (diagnostic) fields
    call create_gungho_prognostics( mesh,                             &
                                    model_data%depository,            &
                                    model_data%prognostic_fields,     &
                                    model_data%diagnostic_fields,     &
                                    model_data%adv_fields_all_outer,  &
                                    model_data%adv_fields_last_outer, &
                                    model_data%mr,                    &
                                    model_data%moist_dyn )

    if (limited_area) call create_lbc_fields( mesh,                         &
                                              model_data%depository,        &
                                              model_data%prognostic_fields, &
                                              model_data%lbc_fields,        &
                                              model_data%lbc_times_list )

    ! Create prognostics used by physics
    if (use_physics) then
      call create_physics_prognostics( mesh, twod_mesh,                  &
                                       model_clock,                      &
                                       model_data%depository,            &
                                       model_data%prognostic_fields,     &
                                       model_data%adv_fields_all_outer,  &
                                       model_data%adv_fields_last_outer, &
                                       model_data%derived_fields,        &
                                       model_data%radiation_fields,      &
                                       model_data%microphysics_fields,   &
                                       model_data%electric_fields,       &
                                       model_data%orography_fields,      &
                                       model_data%turbulence_fields,     &
                                       model_data%convection_fields,     &
                                       model_data%cloud_fields,          &
                                       model_data%surface_fields,        &
                                       model_data%soil_fields,           &
                                       model_data%snow_fields,           &
                                       model_data%chemistry_fields,      &
                                       model_data%aerosol_fields,        &
                                       model_data%stph_fields  )

#ifdef UM_PHYSICS
      ! Create FD prognostic fields
      select case ( prognostic_init_choice )
        case ( init_option_fd_start_dump )
          call create_fd_prognostics(mesh, twod_mesh,       &
                                     model_data%fd_fields,  &
                                     model_data%depository)
      end select

      ! Create and populate collection of fields to be read from ancillary files
      select case ( ancil_choice )
        case ( ancil_option_fixed, ancil_option_updating )
          call create_fd_ancils( model_data%depository,   &
                                 model_data%ancil_fields, &
                                 mesh, twod_mesh ,        &
                                 aerosol_mesh,            &
                                 aerosol_twod_mesh,       &
                                 model_data%ancil_times_list )
      end select
#endif
    end if

  end subroutine create_model_data

  !-------------------------------------------------------------------------------
  !> @brief Initialises the working data set dependent of namelist configuration
  !> @param [inout] model_data The working data set for a model run
  subroutine initialise_model_data( model_data, model_clock, mesh, twod_mesh )

    implicit none

    type( model_data_type ), intent(inout) :: model_data
    class(model_clock_type), intent(in)    :: model_clock
    type( mesh_type ), intent(in), pointer :: mesh
    type( mesh_type ), intent(in), pointer :: twod_mesh

    type( field_type ), pointer :: theta => null()
    type( field_type ), pointer :: rho   => null()
    type( field_type ), pointer :: u     => null()
    type( field_type ), pointer :: exner => null()
    type( field_type ), pointer :: temp_correction_field => null()
    type( field_type ), pointer :: accumulated_fluxes    => null()

    ! Initialise all the physics fields here. We'll then re initialise
    ! them below if need be
    if (use_physics) then
          call init_physics_prognostics_alg( model_data%radiation_fields,    &
                                             model_data%microphysics_fields, &
                                             model_data%electric_fields,     &
                                             model_data%orography_fields,    &
                                             model_data%turbulence_fields,   &
                                             model_data%convection_fields,   &
                                             model_data%cloud_fields,        &
                                             model_data%surface_fields,      &
                                             model_data%soil_fields,         &
                                             model_data%snow_fields,         &
                                             model_data%chemistry_fields,    &
                                             model_data%aerosol_fields,      &
                                             model_data%stph_fields)

    end if

    ! Initialise prognostic fields appropriately
    select case ( prognostic_init_choice )

      case ( init_option_analytic )

        ! Initialise prognostics analytically according to
        ! namelist options

        call init_gungho_prognostics( model_data%prognostic_fields, &
                                      model_data%mr,                &
                                      model_data%moist_dyn,         &
                                      model_data%adv_fields_last_outer )

      case ( init_option_checkpoint_dump )

        ! Initialize prognostics using a checkpoint file
        ! from a previous run

        call read_checkpoint( model_data%prognostic_fields, &
                              model_clock%get_first_step() - 1,   &
                              checkpoint_stem_name )

        ! Update factors for moist dynamics
        call moist_dyn_factors_alg( model_data%moist_dyn, model_data%mr )

      case ( init_option_fd_start_dump )

        if (use_physics) then

          ! Initialise FD prognostic fields from a UM2LFRic dump

          ! Read in from a UM2LFRic dump file
          call init_fd_prognostics_dump( model_data%fd_fields )

          ! Populate prognostics from input finite difference fields
          call map_fd_to_prognostics( model_data%prognostic_fields,          &
                                      model_data%mr,                         &
                                      model_data%moist_dyn,                  &
                                      model_data%fd_fields )

        else
          call log_event("Gungho: Prognostic initialisation from an FD dump not valid "// &
                          "if use_physics is .false., stopping program! ",LOG_LEVEL_ERROR)

        end if

      case ( init_option_fe_start_dump )
        ! Initialise FE prognostic fields from an FE dump
        ! Not yet supported
        call log_event("Gungho: Prognostic initialisation from an FE dump not yet supported, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
      case default
        ! No valid initialisation option selected
        call log_event("Gungho: No valid prognostic initialisation option selected, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)

    end select

    ! initialise coupling fields
#ifdef COUPLED
    if(l_esm_couple) call cpl_init_fields(model_data%cpl_rcv)
#endif

    if (limited_area) then

      select case( lbc_option )
        case ( lbc_option_analytic )
          call init_lbcs_analytic_alg( model_data%prognostic_fields, &
                                       model_data%lbc_fields )

        case ( lbc_option_gungho_file )
          call init_lbcs_file_alg( model_data%lbc_times_list, &
                                   model_clock,                     &
                                   model_data%lbc_fields )

        case ( lbc_option_um2lfric_file )
          call init_lbcs_file_alg( model_data%lbc_times_list, &
                                   model_clock,                     &
                                   model_data%lbc_fields )

        case ( lbc_option_none )
          call log_event( "Gungho: No LBC option specified, yet limited area", LOG_LEVEL_ERROR )
      end select

    end if

    ! Set a default value for temperature_correction_rate for all runs
    model_data%temperature_correction_rate = 0.0_r_def

    ! Assuming this is only relevant for physics runs at the moment
    if (use_physics) then

      ! Initial derived fields
      call model_data%prognostic_fields%get_field('u',u)
      call model_data%prognostic_fields%get_field('theta',theta)
      call model_data%prognostic_fields%get_field('exner',exner)
      call model_data%prognostic_fields%get_field('rho',rho)

      if ( encorr_usage /= encorr_usage_none ) then
        if ( geometry /= geometry_spherical ) then
          write(log_scratch_space, '(a)') &
          'Energy correction valid for spherical geometry only.'
          call log_event(log_scratch_space, log_level_error)
        end if
        call model_data%derived_fields%get_field('temp_correction_field',temp_correction_field)
        call model_data%derived_fields%get_field('accumulated_fluxes',accumulated_fluxes)
      end if

      if (model_clock%is_spinning_up()) then
        call set_any_dof_alg(u, T, 0.0_r_def)
      end if
      call map_physics_fields_alg(u, exner, rho, theta,     &
                                  model_data%moist_dyn,     &
                                  model_data%derived_fields)

      ! Initialise ancillary fields
      select case ( ancil_choice )
        case ( ancil_option_none )
          call log_event( "Gungho: No ancillaries to be read for this run.", LOG_LEVEL_INFO )
#ifdef UM_PHYSICS
        case ( ancil_option_start_dump )
          call log_event( "Gungho: Ancillaries are being read from start dump ", LOG_LEVEL_INFO )
          ! Update the tiled surface temperature with the calculated tstar
          call update_tstar_alg(model_data%surface_fields, &
                                model_data%fd_fields, put_field=.true. )
        case ( ancil_option_fixed, ancil_option_updating )
          call init_variable_fields( model_data%ancil_times_list, &
                                     model_clock, model_data%ancil_fields )
          if (.not. checkpoint_read) then
            ! These parts only need to happen on a new run
            call log_event( "Gungho: Reading ancillaries from file ", LOG_LEVEL_INFO )
            call read_state( model_data%ancil_fields )
            call process_inputs_alg( model_data%ancil_fields,   &
                                     model_data%fd_fields,      &
                                     model_data%surface_fields, &
                                     model_data%soil_fields,    &
                                     model_data%snow_fields,    &
                                     model_data%radiation_fields )

            ! Free up any prognostics not required
            call model_data%depository%remove_field("can_water_in")
            call model_data%fd_fields%remove_field("can_water_in")
            call model_data%depository%remove_field("land_tile_temp")
            call model_data%fd_fields%remove_field("land_tile_temp")
            call model_data%depository%remove_field("tstar_sea_ice")
            call model_data%fd_fields%remove_field("tstar_sea_ice")
            call model_data%depository%remove_field("tile_snow_mass_in")
            call model_data%fd_fields%remove_field("tile_snow_mass_in")
            call model_data%depository%remove_field("n_snow_layers_in")
            call model_data%fd_fields%remove_field("n_snow_layers_in")
            call model_data%depository%remove_field("snow_depth_in")
            call model_data%fd_fields%remove_field("snow_depth_in")
            call model_data%depository%remove_field("tile_snow_rgrain_in")
            call model_data%fd_fields%remove_field("tile_snow_rgrain_in")
            call model_data%depository%remove_field("snow_under_canopy_in")
            call model_data%fd_fields%remove_field("snow_under_canopy_in")
            call model_data%depository%remove_field("snowpack_density_in")
            call model_data%fd_fields%remove_field("snowpack_density_in")

          end if

          ! Free up any ancils which are no longer needed
          call model_data%depository%remove_field("land_area_fraction")
          call model_data%ancil_fields%remove_field("land_area_fraction")
          call model_data%depository%remove_field("land_tile_fraction")
          call model_data%ancil_fields%remove_field("land_tile_fraction")
          call model_data%depository%remove_field("stdev_topog_index")
          call model_data%ancil_fields%remove_field("stdev_topog_index")
#endif
        case default
          ! No valid ancil option selected
          call log_event("Gungho: No valid ancillary initialisation option selected, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
      end select

      ! Initialise energy correction
      if ( encorr_usage /= encorr_usage_none ) then
        call compute_total_mass_alg( model_data%total_dry_mass, rho, mesh )

        call compute_total_energy_alg( model_data%total_energy_previous,                &
                                       model_data%derived_fields, u, theta, exner, rho, &
                                       model_data%mr,                                   &
                                       mesh, twod_mesh )

        ! Initialise flux sum to zero
        call scalar_to_field_alg(0.0_r_def, accumulated_fluxes)

        if ( checkpoint_read ) then
          ! Read scalar temperature_correction_rate from checkponted field
          ! temp_correction_field
          call field_to_scalar_alg(model_data%temperature_correction_rate, &
                                   temp_correction_field)
          write(log_scratch_space, '(''restart_temperature_correction_rate = '', e30.22)') &
               model_data%temperature_correction_rate
          call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
        else
          model_data%temperature_correction_rate = 0.0_r_def
          call scalar_to_field_alg(0.0_r_def, temp_correction_field)
        end if

      end if

    end if

  end subroutine initialise_model_data

  !----------------------------------------------------------------------------
  !> @brief Writes out a checkpoint and dump file dependent on namelist
  !> options
  !> @param[inout] model_data The working data set for the model run
  !> @param[in] model_clock Model time.
  !>
  subroutine output_model_data( model_data, model_clock )

    implicit none

    type( model_data_type ), intent(inout), target :: model_data
    class(model_clock_type), intent(in)            :: model_clock

    type( field_collection_type ), pointer :: fd_fields => null()
    type( field_collection_type ), pointer :: prognostic_fields => null()

    ! Get pointers to field collections for use downstream
    fd_fields => model_data%fd_fields
    prognostic_fields => model_data%prognostic_fields

    !===================== Write fields to dump ======================!
    if( use_physics ) then

      ! Current dump writing is only relevant for physics runs at the moment
      if (write_dump) then

        call log_event("Gungho: writing FD fields to dump is not yet supported", LOG_LEVEL_ERROR)

        ! Write prognostic fields to dump
        call write_state(fd_fields)

      end if

    end if

    !=================== Write fields to checkpoint files ====================!
    if( checkpoint_write ) then
       call write_checkpoint( prognostic_fields, model_clock, checkpoint_stem_name )
    end if

  end subroutine output_model_data

  !>@brief Routine to destroy all the field collections in the working data set
  !> @param[inout] model_data The working data set for a model run
  subroutine finalise_model_data( model_data )

    implicit none

      type(model_data_type), intent(inout) :: model_data

      ! Clear all the fields in each field collection
      call model_data%depository%clear()
      call model_data%prognostic_fields%clear()
      call model_data%diagnostic_fields%clear()
      call model_data%adv_fields_last_outer%clear()
      call model_data%adv_fields_all_outer%clear()
      call model_data%derived_fields%clear()
      call model_data%radiation_fields%clear()
      call model_data%microphysics_fields%clear()
      call model_data%electric_fields%clear()
      call model_data%orography_fields%clear()
      call model_data%turbulence_fields%clear()
      call model_data%convection_fields%clear()
      call model_data%cloud_fields%clear()
      call model_data%surface_fields%clear()
      call model_data%soil_fields%clear()
      call model_data%snow_fields%clear()
      call model_data%chemistry_fields%clear()
      call model_data%aerosol_fields%clear()
      call model_data%stph_fields%clear()
      call model_data%fd_fields%clear()
#ifdef COUPLED
      if(l_esm_couple) then
         call model_data%cpl_snd%clear()
         call model_data%cpl_rcv%clear()
      endif
#endif
      call model_data%lbc_fields%clear()
      call model_data%ls_fields%clear()
      if (allocated(model_data%mr)) deallocate(model_data%mr)
      if (allocated(model_data%moist_dyn)) deallocate(model_data%moist_dyn)
      if (allocated(model_data%ls_mr)) deallocate(model_data%ls_mr)
      if (allocated(model_data%ls_moist_dyn)) deallocate(model_data%ls_moist_dyn)

      call log_event( 'finalise_model_data: all fields have been cleared', &
                       LOG_LEVEL_INFO )

  end subroutine finalise_model_data

end module gungho_model_data_mod
