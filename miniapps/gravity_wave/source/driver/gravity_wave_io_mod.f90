!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Controls output (diags/checkpointing) related information used by
!>        the model

module gravity_wave_io_mod

  use clock_mod,               only : clock_type
  use constants_mod,           only : i_def, i_native
  use field_mod,               only : field_type
  use io_config_mod,           only : use_xios_io
  use lfric_xios_io_mod,       only : initialise_xios
  use xios,                    only : xios_context_finalize, &
                                      xios_update_calendar

  implicit none

  private
  public initialise_io, finalise_io

contains

  !> @brief Initialises output (diags/checkpointing) used by the model
  !> @param [inout] comm The MPI communicator for use within the model
  !> @param [in] clock Model time.
  !> @param [in] mesh_id The identifier of the primary mesh
  !> @param [in] twod_mesh_id The identifier of the primary 2d mesh
  !> @param [in] chi_xyz A size 3 array of fields holding the (X,Y,Z)
  !>                     coordinates of the mesh.
  !> @param [in] xios_ctx XIOS context identifier
  !>
  subroutine initialise_io(comm, clock, mesh_id, twod_mesh_id, chi_xyz, xios_ctx)

    implicit none

    integer(i_native), intent(in) :: comm
    type(clock_type),  intent(in) :: clock
    integer(i_def),    intent(in) :: mesh_id, twod_mesh_id
    type(field_type),  intent(in) :: chi_xyz(3)
    character(len=*),  intent(in) :: xios_ctx

  !----------------------------------------------------------------------------
  ! IO init
  !----------------------------------------------------------------------------

  ! If using XIOS for diagnostic output or checkpointing, then set up XIOS
  ! domain and context
  if ( use_xios_io ) then

    call initialise_xios( xios_ctx,     &
                          comm,         &
                          clock,        &
                          mesh_id,      &
                          twod_mesh_id, &
                          chi_xyz )

    if (clock%is_initialisation()) then
      ! Make sure XIOS calendar is set to timestep 1 as it starts there
      ! not timestep 0.
      call xios_update_calendar(1)
    end if

  end if

  end subroutine initialise_io

  !> @brief Finalises output related functions used by the model
  subroutine finalise_io()

    implicit none

    ! Finalise XIOS context if we used it for diagnostic output or checkpointing
    if ( use_xios_io ) then
      call xios_context_finalize()
    end if

  end subroutine finalise_io

end module gravity_wave_io_mod
