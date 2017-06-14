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
  use constants_mod

  implicit none

  ! The DEPENDS ON statment below results in a dependency on the 
  ! UM code and so the build system should build in the
  ! relevant dependencies
  !DEPENDS ON: qsat

  character(len=100) :: fname = 'umphysics_testbuild-checksums.txt'
  integer            :: stat

  open( 9, file=fname, status="replace", iostat=stat)
  if (stat /= 0) then
    print*, "Unable to open checksum file"
  end if
  write(9, '(A)' ) 'UM physics build test program successful'

end program umphysics_testbuild
