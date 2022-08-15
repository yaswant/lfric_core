!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> Da_Dev miniapp program support functions.
!>
!> Originally these were "block" constructs within the program but neither
!> GNU or Intel Fortran where properly able to cope with that.
!>
module da_dev_mod

  use log_mod, only : log_event,         &
                      log_scratch_space, &
                      LOG_LEVEL_ALWAYS,  &
                      LOG_LEVEL_ERROR


  implicit none

  private
  public :: load_configuration

  contains

  !> Loads run-time configuration and ensures everything is ship-shape.
  !>
  subroutine load_configuration( filename, model_name )

    use configuration_mod, only : read_configuration, &
                                  ensure_configuration

    implicit none

    character(*), intent(in) :: filename
    character(*), intent(in) :: model_name

    character(*), parameter ::                           &
        required_configuration(5) =  [ 'base_mesh     ', &
                                       'extrusion     ', &
                                       'finite_element', &
                                       'partitioning  ', &
                                       'planet        ']

    logical              :: okay
    logical, allocatable :: success_map(:)
    integer              :: i

    allocate( success_map(size(required_configuration)) )

    call log_event( 'Loading '//model_name//' configuration ...', &
                    LOG_LEVEL_ALWAYS )

    call read_configuration( filename )

    okay = ensure_configuration( required_configuration, success_map )
    if (.not. okay) then
      write( log_scratch_space, '(A)' ) &
                             'The following required namelists were not loaded:'
      do i = 1,size(required_configuration)
        if (.not. success_map(i)) &
          log_scratch_space = trim(log_scratch_space) // ' ' &
                              // required_configuration(i)
      end do
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    deallocate( success_map )

  end subroutine load_configuration

end module da_dev_mod
