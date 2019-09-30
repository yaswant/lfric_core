!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Steps the gungho app through one timestep

!> @details Handles the stepping (for a single timestep) of the
!>          gungho app

module gungho_step_mod

  use conservation_algorithm_mod,     only : conservation_algorithm
  use constants_mod,                  only : i_def
  use diagnostics_calc_mod,           only : write_density_diagnostic
  use field_collection_mod,           only : field_collection_type
  use field_mod,                      only : field_type
  use formulation_config_mod,         only : transport_only
  use io_config_mod,                  only : write_diag, &
                                             write_minmax_tseries
  use log_mod,                        only : LOG_LEVEL_INFO
  use minmax_tseries_mod,             only : minmax_tseries
  use mr_indices_mod,                 only : nummr
  use iter_timestep_alg_mod,          only : iter_alg_step
  use moist_dyn_mod,                  only : num_moist_factors
  use rk_alg_timestep_mod,            only : rk_alg_step
  use rk_transport_mod,               only : rk_transport_step
  use timestepping_config_mod,        only : method, &
                                             method_semi_implicit, &
                                             method_rk
  use transport_config_mod,           only : scheme, &
                                             scheme_method_of_lines

  implicit none

  private
  public step

  contains

  !> @brief Steps the gungho app through one timestep
  !> @param[in] mesh_id The identifier of the primary mesh
  !> @param[in] twod_mesh_id The identifier of the two-dimensional mesh
  !> @param[inout] prognostic_fields A collection of all the prognostic fields
  !> @param[inout] diagnostic_fields A collection of all the diagnostic fields
  !> @param[inout] mr Array of fields containing the mixing ratios
  !> @param[inout] moist_dyn Array of fields containing factors for moist dynamics
  !> @param[inout] derived_fields Collection of FD fields derived from FE fields
  !> @param[inout] cloud_fields Collection of cloud fields
  !> @param[inout] twod_fields 2D field collection for physics
  !> @param[inout] radstep_fields Collection of radiation timestep fields
  !> @param[inout] physics_incs Collection of physics increments
  !> @param[inout] jules_ancils Collection of Jules ancillaries
  !> @param[inout] jules_prognostics Collection of Jules prognostics
  !> @param[in] timestep number of current timestep
  subroutine step(mesh_id,           &
                  twod_mesh_id,      &
                  prognostic_fields, &
                  diagnostic_fields, &
                  mr,                &
                  moist_dyn,         &
                  derived_fields,    &
                  cloud_fields,      &
                  twod_fields,       &
                  radstep_fields,    &
                  physics_incs,      &
                  jules_ancils,      &
                  jules_prognostics, &
                  timestep)

    implicit none

    integer(i_def),                intent(in)    :: mesh_id
    integer(i_def),                intent(in)    :: twod_mesh_id
    type( field_collection_type ), intent(inout) :: prognostic_fields
    type( field_collection_type ), intent(inout) :: diagnostic_fields
    type( field_type ),            intent(inout) :: mr(nummr)
    type( field_type ),            intent(inout) :: moist_dyn(num_moist_factors)
    type( field_collection_type ), intent(inout) :: derived_fields
    type( field_collection_type ), intent(inout) :: cloud_fields
    type( field_collection_type ), intent(inout) :: twod_fields
    type( field_collection_type ), intent(inout) :: radstep_fields
    type( field_collection_type ), intent(inout) :: physics_incs
    type( field_collection_type ), intent(inout) :: jules_ancils
    type( field_collection_type ), intent(inout) :: jules_prognostics
    integer(i_def),                intent(in)    :: timestep

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()

    ! Get pointers to fields in the prognostic/diagnostic field collections
    ! for use downstream
    theta => prognostic_fields%get_field('theta')
    u => prognostic_fields%get_field('u')
    rho => prognostic_fields%get_field('rho')
    exner => prognostic_fields%get_field('exner')

    if ( transport_only ) then
      select case( scheme )
        case ( scheme_method_of_lines )
          call rk_transport_step( u, rho, theta)
      end select
      call write_density_diagnostic(rho, timestep)
      if ( write_diag ) &
        call conservation_algorithm(timestep, rho, u, theta, exner)
    else  ! Not transport_only
      select case( method )
        case( method_semi_implicit )  ! Semi-Implicit
          call iter_alg_step(u, rho, theta, exner, mr, moist_dyn,           &
                             derived_fields, cloud_fields, twod_fields,     &
                             radstep_fields, physics_incs,                  &
                             jules_ancils, jules_prognostics,               &
                             timestep, twod_mesh_id)
        case( method_rk )             ! RK
          call rk_alg_step(u, rho, theta, moist_dyn, exner)
      end select

      if ( write_diag ) call conservation_algorithm(timestep, rho, u, theta, exner)

      if(write_minmax_tseries) call minmax_tseries(u, 'u', mesh_id)

      call u%log_minmax(LOG_LEVEL_INFO, ' u')

    end if

  end subroutine step

end module gungho_step_mod
