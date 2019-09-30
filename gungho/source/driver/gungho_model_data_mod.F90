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
  use log_mod,                    only : log_event, &
                                         LOG_LEVEL_INFO

  implicit none

  private

  !> Holds the working data set for a model run.
  !>
  type :: model_data_type

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

  end type model_data_type

  public model_data_type, finalise_model_data

contains

  !>@brief Routine to destroy all the field collections in the working data set
  !> @param[inout] model_data The working data set for a model run
  subroutine finalise_model_data( model_data )
 
    implicit none

      type(model_data_type), intent(inout) :: model_data

      ! Clear all the fields in each field collection
      call model_data%depository%clear()
      call model_data%prognostic_fields%clear()
      call model_data%diagnostic_fields%clear()
      call model_data%jules_ancils%clear()
      call model_data%jules_prognostics%clear()
      call model_data%derived_fields%clear()
      call model_data%cloud_fields%clear()
      call model_data%twod_fields%clear()
      call model_data%radstep_fields%clear()
      call model_data%physics_incs%clear()
      call model_data%fd_fields%clear()

      call log_event( 'finalise_model_data: all fields have been cleared', LOG_LEVEL_INFO )

  end subroutine finalise_model_data

end module gungho_model_data_mod
