!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> Set up and tear down for whole-suite fixtures.
!>
module project_suite_fixture_mod

  use ESMF
  use constants_mod, only: str_def, i_def

  implicit none
  private
  public :: get_esmf_handle, set_up_suite, tear_down_suite

  type(ESMF_VM) :: vm

contains

  !> Get the ESMF virtual machine handle.
  !>
  type(ESMF_VM) function get_esmf_handle()

    implicit none

    get_esmf_handle = vm

  end function get_esmf_handle

  !> Set up any whole-suite fixtures.
  !> This is called by pFUnit, not developer code.
  !>
  subroutine set_up_suite()

    implicit none

    integer(i_def)     :: rc
    character(str_def) :: project_name

    CALL get_environment_variable("PROJECT_NAME", project_name, status=rc)

    if (rc <= 0) then
      project_name = trim(project_name)//'-'
    else
      project_name = ''
    end if

    call ESMF_Initialize( vm=vm, &
                          defaultlogfilename=trim(project_name)//"pfunit.log", &
                          logkindflag=ESMF_LOGKIND_MULTI, &
                          rc=rc )
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize( endflag=ESMF_END_ABORT )

  end subroutine set_up_suite

  !> Tear down any whole-suite fixtures.
  !> This is called by pFUnit, not developer code.
  !>
  subroutine tear_down_suite()

    implicit none

    integer(i_def) :: rc

    ! MPI must be kept open so that subsequent finalisation can close it.
    !
    call ESMF_Finalize( endflag=ESMF_END_KEEPMPI, rc=rc )

  end subroutine tear_down_suite

end module project_suite_fixture_mod
