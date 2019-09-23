!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Container for Gung Ho model run working data set.
module gungho_model_data_mod

  use mr_indices_mod,             only : nummr
  use moist_dyn_mod,              only : num_moist_factors
  use field_mod,                  only : field_type
  use field_collection_mod,       only : field_collection_type

  implicit none

  private

  !> Holds the working data set for a model run.
  !>
  type, public :: gungho_model_data_type

    private

    !> Stores all the fields used by the model - the remaining collections store
    !> pointers to the data in the depository.
    type( field_collection_type ), public :: depository

    !> @name Fields needed to time-step the model.
    !> @{

    !> Ancillary fields for Jules
    type( field_collection_type ), public :: jules_ancils
    !> Prognostic fields for Jules
    type( field_collection_type ), public :: jules_prognostics
    !> All the prognostic fields (except for field arrays: auxiliary prognostic)
    type( field_collection_type ), public :: prognostic_fields
    !> All the diagnostic fields
    type( field_collection_type ), public :: diagnostic_fields
    !> FD fields derived from FE fields for use in physics time-stepping schemes
    type( field_collection_type ), public :: derived_fields
    !> Cloud fields used by the physics time-stepping schemes
    type( field_collection_type ), public :: cloud_fields
    !> 2D fields used by the UM physics
    type( field_collection_type ), public :: twod_fields
    !> Fields used for the radiation timestep
    type( field_collection_type ), public :: radstep_fields
    !> Increment fields used by the physics time-stepping schemes
    type( field_collection_type ), public :: physics_incs
    !> Array of fields containing the moisture mixing ratios (auxiliary prognostic)
    type( field_type ), public :: mr(nummr)
    !> Array of fields containing the moist dynamics (auxiliary prognostic)
    type( field_type ), public :: moist_dyn(num_moist_factors)
    !> @}

    !> FD fields used to read initial conditions from LFRic-Input files
    type( field_collection_type ), public :: fd_fields

  end type gungho_model_data_type

end module gungho_model_data_mod
