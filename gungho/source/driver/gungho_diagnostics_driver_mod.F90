!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Outputs diagnostics from gungho/lfric_atm

!> @details Calls the routine that generates diagnostic output for
!>          gungho/lfric_atm. This is only a temporary
!>          hard-coded solution in lieu of a proper dianostic system

module gungho_diagnostics_driver_mod

  use clock_mod,                 only : clock_type
  use constants_mod,             only : i_def, str_def
  use boundaries_config_mod,     only : limited_area, output_lbcs
  use diagnostics_io_mod,        only : write_scalar_diagnostic, &
                                        write_vector_diagnostic
  use diagnostics_calc_mod,      only : write_divergence_diagnostic, &
                                        write_hydbal_diagnostic, &
                                        write_vorticity_diagnostic
  use field_collection_iterator_mod, &
                                 only : field_collection_iterator_type
  use field_collection_mod,      only : field_collection_type
  use diagnostic_alg_mod,        only : column_total_diagnostics_alg, &
                                        calc_wbig_diagnostic_alg
  use gungho_model_data_mod,     only : model_data_type
  use field_mod,                 only : field_type
  use field_parent_mod,          only : field_parent_type
  use formulation_config_mod,    only : use_moisture, &
                                        use_physics
  use fs_continuity_mod,         only : W3, Wtheta
  use integer_field_mod,         only : integer_field_type
  use moist_dyn_mod,             only : num_moist_factors
  use mr_indices_mod,            only : nummr, mr_names
  use section_choice_config_mod, only : cloud, cloud_um
  use log_mod,                   only : log_event, &
                                        LOG_LEVEL_INFO
  use geometric_constants_mod,   only : get_panel_id

  implicit none

  private
  public gungho_diagnostics_driver

contains

  !> @brief Outputs simple diagnostics from Gungho/LFRic
  !> @param[in] mesh_id The identifier of the primary mesh
  !> @param[in] twod_mesh_id The identifier given to the current 2d mesh
  !> @param[in] model_data The working data set for the model run
  !> @param[in] timestep The timestep at which the fields are valid
  !> @param[in] nodal_output_on_w3 Flag that determines if vector fields
  !>                  should be projected to W3 for nodal output
  subroutine gungho_diagnostics_driver( mesh_id,      &
                                        twod_mesh_id, &
                                        model_data,   &
                                        clock,        &
                                        nodal_output_on_w3 )

    implicit none

    integer(i_def),        intent(in)         :: mesh_id, twod_mesh_id
    type(model_data_type), intent(in), target :: model_data
    class(clock_type),     intent(in)         :: clock
    logical,               intent(in)         :: nodal_output_on_w3

    type( field_collection_type ), pointer :: prognostic_fields => null()
    type( field_collection_type ), pointer :: diagnostic_fields => null()
    type( field_collection_type ), pointer :: lbc_fields => null()
    type( field_type ),            pointer :: mr(:) => null()
    type( field_type ),            pointer :: moist_dyn(:) => null()
    type( field_collection_type ), pointer :: derived_fields => null()
    type( field_collection_type ), pointer :: cloud_fields => null()

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: h_u => null()
    type( field_type), pointer :: v_u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()
    type( field_type), pointer :: panel_id => null()
    type( field_type), pointer :: u_star => null()
    type( field_type), pointer :: w_physics => null()
    ! Iterator for field collection
    type(field_collection_iterator_type)  :: iterator

    ! A pointer used for retrieving fields from collections
    ! when iterating over them
    class( field_parent_type ), pointer :: field_ptr  => null()

    character(str_def) :: name

    integer :: i, fs

    call log_event("Gungho: writing diagnostic output", LOG_LEVEL_INFO)

    ! Get pointers to field collections for use downstream
    prognostic_fields => model_data%prognostic_fields
    diagnostic_fields => model_data%diagnostic_fields
    lbc_fields => model_data%lbc_fields
    mr => model_data%mr
    moist_dyn => model_data%moist_dyn
    derived_fields => model_data%derived_fields
    cloud_fields => model_data%cloud_fields
    panel_id => get_panel_id(mesh_id)

    ! Can't just iterate through the prognostic/diagnostic collections as
    ! some fields are scalars and some fields are vectors, so explicitly
    ! extract all fields from the collections and output each of them in turn
    theta => prognostic_fields%get_field('theta')
    u => prognostic_fields%get_field('u')
    rho => prognostic_fields%get_field('rho')
    exner => prognostic_fields%get_field('exner')

    ! Scalar fields
    call write_scalar_diagnostic('rho', rho, &
                                 clock, mesh_id, nodal_output_on_w3)
    call write_scalar_diagnostic('theta', theta, &
                                 clock, mesh_id, nodal_output_on_w3)
    call write_scalar_diagnostic('exner', exner, &
                                 clock, mesh_id, nodal_output_on_w3)

    ! Vector fields
    call write_vector_diagnostic('u', u, &
                                 clock, mesh_id, nodal_output_on_w3)
    call write_vorticity_diagnostic( u, clock )

    ! Moisture fields
    if (use_moisture) then
      do i=1,nummr
        call write_scalar_diagnostic( trim(mr_names(i)), mr(i), &
                                      clock, mesh_id, nodal_output_on_w3 )
      end do
    end if

    if (limited_area) then
      if (output_lbcs) then
        theta => lbc_fields%get_field('lbc_theta')
        u => lbc_fields%get_field('lbc_u')
        rho => lbc_fields%get_field('lbc_rho')
        exner => lbc_fields%get_field('lbc_exner')

        h_u => lbc_fields%get_field('lbc_h_u')
        v_u => lbc_fields%get_field('lbc_v_u')

        ! Scalar fields
        call write_scalar_diagnostic('lbc_rho', rho, &
                                 clock, mesh_id, nodal_output_on_w3)
        call write_scalar_diagnostic('lbc_theta', theta, &
                                 clock, mesh_id, nodal_output_on_w3)
        call write_scalar_diagnostic('lbc_exner', exner, &
                                 clock, mesh_id, nodal_output_on_w3)
        call write_scalar_diagnostic('readlbc_v_u', v_u, &
                                 clock, mesh_id, nodal_output_on_w3)

        ! Vector fields
        call write_vector_diagnostic('lbc_u', u, &
                                 clock, mesh_id, nodal_output_on_w3)
        call write_vector_diagnostic('readlbc_h_u', h_u, &
                                 clock, mesh_id, nodal_output_on_w3)
      endif
    endif

    ! Cloud fields - are all scalars - so iterate over collection
    if (use_physics .and. cloud == cloud_um) then

      call iterator%initialise(cloud_fields)
      do
        if ( .not.iterator%has_next() ) exit
        field_ptr => iterator%next()

        select type(field_ptr)
          type is (field_type)
            name = trim(adjustl( field_ptr%get_name() ))
            call write_scalar_diagnostic( trim(name), field_ptr, &
                                          clock,                 &
                                          mesh_id, nodal_output_on_w3 )
        end select
      end do
      field_ptr => null()
    end if

    ! Derived physics fields (only those on W3 or Wtheta)
    if (use_physics .and. .not. clock%is_initialisation()) then

      call iterator%initialise(derived_fields)
      do
        if ( .not.iterator%has_next() ) exit
        field_ptr => iterator%next()
        select type(field_ptr)
          type is (field_type)
            fs = field_ptr%which_function_space()
            if ( fs == W3 .or. fs == Wtheta ) then
              name = trim(adjustl( field_ptr%get_name() ))
              call write_scalar_diagnostic( trim(name), field_ptr, &
                                            clock,                 &
                                            mesh_id, nodal_output_on_w3 )
            end if
        end select
      end do
      field_ptr => null()

      ! Output u_star as diagnostic
      u_star => derived_fields%get_field('u_star')
      call write_vector_diagnostic('u_star',u_star,clock,mesh_id,nodal_output_on_w3)

      ! Get w_physics for WBig calculation
      w_physics => derived_fields%get_field('w_physics')
      call calc_wbig_diagnostic_alg(w_physics, mesh_id)

    end if

    ! Other derived diagnostics with special pre-processing
    call write_divergence_diagnostic(u, clock, mesh_id)
    call write_hydbal_diagnostic(theta, moist_dyn, exner, mesh_id)
    call column_total_diagnostics_alg(rho, mr, mesh_id, twod_mesh_id)


  end subroutine gungho_diagnostics_driver

end module gungho_diagnostics_driver_mod
