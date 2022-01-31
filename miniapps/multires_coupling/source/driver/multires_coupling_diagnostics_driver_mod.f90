!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Outputs the diagnostics from the multires_coupling miniapp

!> @details Calls the routine that generates diagnostic output for all
!>          fields used by the multires_coupling miniapp.

module multires_coupling_diagnostics_driver_mod

  use clock_mod,                          only : clock_type
  use constants_mod,                      only : i_def
  use field_mod,                          only : field_type
  use field_collection_mod,               only : field_collection_type
  use formulation_config_mod,             only : use_physics
  use gungho_model_data_mod,              only : model_data_type
  use gungho_diagnostics_driver_mod,      only : gungho_diagnostics_driver
  use map_physics_fields_alg_mod,         only : map_physics_fields_alg

  implicit none

  private
  public multires_coupling_diagnostics_driver

contains

  !> @brief Outputs the diagnostics from the multires_coupling miniapp.
  !> @param [in] dynamics_mesh_id The identifier of the dynamics mesh
  !> @param [in] dynamics_2D_mesh_id The identifier of the dynamics 2D mesh
  !> @param [in,out] dynamics_model_data A collection containing the fields on the
  !>                dynamics mesh
  !> @param [in] clock Model time.
  !> @param [in] W3_project Flag that determines if vector fields should be
  !>                        projected to W3
  subroutine multires_coupling_diagnostics_driver( dynamics_mesh_id,       &
                                                   dynamics_2D_mesh_id,    &
                                                   dynamics_model_data,    &
                                                   clock, W3_project )

    implicit none

    type(model_data_type), intent(inout), target   :: dynamics_model_data
    integer(kind=i_def),   intent(in)              :: dynamics_mesh_id
    integer(kind=i_def),   intent(in)              :: dynamics_2D_mesh_id
    class(clock_type),     intent(in)              :: clock
    logical,               intent(in)              :: W3_project

    type( field_collection_type ), pointer :: dynamics_prognostic_fields => null()
    type( field_collection_type ), pointer :: dynamics_derived_fields => null()
    type( field_type ),            pointer :: dynamics_mr(:) => null()
    type( field_type ),            pointer :: dynamics_moist_dyn(:) => null()

    type(field_type), pointer :: dynamics_u => null()
    type(field_type), pointer :: dynamics_rho => null()
    type(field_type), pointer :: dynamics_theta => null()
    type(field_type), pointer :: dynamics_exner => null()

    dynamics_prognostic_fields => dynamics_model_data%prognostic_fields
    dynamics_derived_fields => dynamics_model_data%derived_fields
    dynamics_mr => dynamics_model_data%mr
    dynamics_moist_dyn => dynamics_model_data%moist_dyn

    ! Can't just iterate through the collection as some fields are scalars
    ! and some fields are vectors, so explicitly extract all fields from
    ! the collection and output each of them
    dynamics_u => dynamics_prognostic_fields%get_field('u')
    dynamics_rho => dynamics_prognostic_fields%get_field('rho')
    dynamics_theta => dynamics_prognostic_fields%get_field('theta')
    dynamics_exner => dynamics_prognostic_fields%get_field('exner')

    if (use_physics) then
      call map_physics_fields_alg(dynamics_u, dynamics_exner,      &
                                  dynamics_rho, dynamics_theta,    &
                                  dynamics_moist_dyn, dynamics_derived_fields)
    end if

    call gungho_diagnostics_driver( dynamics_mesh_id, dynamics_2D_mesh_id, &
                                    dynamics_model_data, clock, W3_project )

  end subroutine multires_coupling_diagnostics_driver

end module multires_coupling_diagnostics_driver_mod
