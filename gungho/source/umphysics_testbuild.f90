!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @mainpage umphysics_testbuild
!> Test program for the build of UM physics source code
!> @brief Test program for the build of UM physics source code.
!>        The program simply has a dependency on the a subset of
!>        the UM physics code, which should then build. However, the
!>        program itself just performs a trivial print statement.
program umphysics_testbuild

  ! The build system requires some dependence so we
  ! Add an arbitrary dependence on constants_mod
  use constants_mod, only: r_def, r_single, r_double

  ! This is a UM module, we arbitrarily chose g see what precision it is
  use planet_constants_mod, only: g
  implicit none

  character(len=100)  :: fname = 'umphysics_testbuild-checksums.txt'
  integer             :: stat
  integer, parameter  :: funit = 9
  real(kind=r_def)    :: x1
  real(kind=r_single) :: x2
  real(kind=r_double) :: x3
  real                :: x4

  open( funit, file=fname, status="replace", iostat=stat)
  if (stat /= 0) then
    print*, "Unable to open checksum file"
  end if
  write(funit, '(A)' ) 'UM physics build test program successful'
  write(funit, '(A40,I13)' ) 'LFRic defined precision reals have kind: ', kind(x1)
  write(funit, '(A40,I13)' ) 'LFRic single precision reals have kind:  ', kind(x2)
  write(funit, '(A40,I13)' ) 'LFRic double precision reals have kind:  ', kind(x3)
  write(funit, '(A40,I13)' ) 'LFRic default precision reals have kind: ', kind(x4)
  write(funit, '(A40,I13)' ) 'UM default precision reals have kind:    ', kind(g)

  close(funit)

end program umphysics_testbuild
