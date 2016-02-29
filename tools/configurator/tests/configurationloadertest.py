#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

import unittest
import StringIO

import configurator.configurationloader as loader

###############################################################################
class LoaderTest( unittest.TestCase ):
    def setUp( self ):
        self.maxDiff = None

    ###########################################################################
    def testEmpty( self ):
        expectedSource = '''
!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
! Handles the loading of namelists.
!
module configuration_mod

  use constants_mod, only : i_native, l_def, str_def, str_max_filename
  use log_mod,       only : log_scratch_space, log_event, LOG_LEVEL_ERROR

  implicit none

  private
  public :: read_configuration, ensure_configuration

  interface read_configuration
    procedure read_configuration_file
    procedure read_configuration_unit
  end interface read_configuration

contains

  ! Reads configuration namelists from a file.
  !
  ! [in] filename File holding the namelists.
  !
  ! TODO: Assumes namelist tags come at the start of lines.
  ! TODO: Support "namelist file" namelists which recursively call this
  !       procedure to load other namelist files.
  !
  subroutine read_configuration_file( filename )

    use io_utility_mod, only : open_file, close_file

    implicit none

    character(*), intent(in) :: filename

    integer(i_native) :: unit
    integer(i_native) :: success

    unit = open_file( filename )
    call read_configuration_unit( unit )
    call close_file( unit )

  end subroutine read_configuration_file

  ! Reads configuration namelists from an I/O unit.
  !
  ! [in] unit File holding the namelists.
  !
  subroutine read_configuration_unit( unit )

    use io_utility_mod, only : read_line

    implicit none

    integer(i_native), intent(in) :: unit

    character(str_def) :: buffer
    logical(l_def)     :: continue
    character(str_max_filename) :: filename

    text_line_loop: do

      continue = read_line( unit, buffer )
      if (.not. continue) exit

      if (buffer(1:1) == '&') then
        select case (trim(buffer(2:)))

          case default
            inquire( unit, name=filename )
            write( log_scratch_space, '(A, A, A, A)' ) &
                 "Unrecognised namelist """,           &
                 trim(buffer),                         &
                 """ found in file ",                  &
                 filename
            call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end select
      end if
    end do text_line_loop

  end subroutine read_configuration_unit

  ! Checks that the requested namelists have been loaded.
  !
  ! [in]  names List of namelists.
  ! [out] success_mask Marks corresponding namelists as having failed.
  !
  ! [return] Overall success.
  !
  function ensure_configuration( names, success_mask )

    implicit none

    character(*),             intent(in)  :: names(:)
    logical(l_def), optional, intent(out) :: success_mask(:)
    logical(l_def)                        :: ensure_configuration

    integer(i_native) :: i
    logical           :: configuration_found = .True.

    if (present(success_mask) &
        .and. (size(success_mask, 1) /= size(names, 1))) then
      call log_event( 'Arguments "name" and "success_mask" to function' &
                      // '"ensure_configuration" are different shapes', &
                      LOG_LEVEL_ERROR )
    end if

    ensure_configuration = .True.

    name_loop: do i = 1, size(names)
      select cae(trim( names(i) ))

        case default
          write( log_scratch_space, '(A, A, A)' )          &
               "Tried to ensure unrecognised namelist """, &
               trim(names(i)),                             &
               """ was loaded"
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end select

      ensure_configuration = ensure_configuration .and. configuration_found

      if (present(success_mask)) success_mask(i) = configuration_found

    end do name_loop

  end function ensure_configuration

end module configuration_mod
        '''.strip().replace( '~', ' ' )

        outputFile = StringIO.StringIO()
        uut = loader.ConfigurationLoader()
        uut.writeModule( outputFile )

        self.assertMultiLineEqual( expectedSource, outputFile.getvalue() )

    ###########################################################################
    def testEmpty( self ):
        expectedSource = '''
!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
! Handles the loading of namelists.
!
module configuration_mod

  use constants_mod, only : i_native, l_def, str_def, str_max_filename
  use log_mod,       only : log_scratch_space, log_event, LOG_LEVEL_ERROR

{usage}

  implicit none

  private
  public :: read_configuration, ensure_configuration

  interface read_configuration
    procedure read_configuration_file
    procedure read_configuration_unit
  end interface read_configuration

contains

  ! Reads configuration namelists from a file.
  !
  ! [in] filename File holding the namelists.
  !
  ! TODO: Assumes namelist tags come at the start of lines.
  ! TODO: Support "namelist file" namelists which recursively call this
  !       procedure to load other namelist files.
  !
  subroutine read_configuration_file( filename )

    use io_utility_mod, only : open_file, close_file

    implicit none

    character(*), intent(in) :: filename

    integer(i_native) :: unit

    unit = open_file( filename )
    call read_configuration_unit( unit )
    call close_file( unit )

  end subroutine read_configuration_file

  ! Reads configuration namelists from an I/O unit.
  !
  ! [in] unit File holding the namelists.
  !
  subroutine read_configuration_unit( unit )

    use io_utility_mod, only : read_line

    implicit none

    integer(i_native), intent(in) :: unit

    character(str_def) :: buffer
    logical(l_def)     :: continue
    character(str_max_filename) :: filename

    text_line_loop: do

      continue = read_line( unit, buffer )
      if (.not. continue) exit

      if (buffer(1:1) == '&') then
        select case (trim(buffer(2:)))
{readSection}
          case default
            inquire( unit, name=filename )
            write( log_scratch_space, '(A, A, A, A)' ) &
                 "Unrecognised namelist """,           &
                 trim(buffer),                         &
                 """ found in file ",                  &
                 filename
            call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end select
      end if
    end do text_line_loop

  end subroutine read_configuration_unit

  ! Checks that the requested namelists have been loaded.
  !
  ! [in]  names List of namelists.
  ! [out] success_mask Marks corresponding namelists as having failed.
  !
  ! [return] Overall success.
  !
  function ensure_configuration( names, success_mask )

    implicit none

    character(*),             intent(in)  :: names(:)
    logical(l_def), optional, intent(out) :: success_mask(:)
    logical(l_def)                        :: ensure_configuration

    integer(i_native) :: i
    logical           :: configuration_found = .True.

    if (present(success_mask) &
        .and. (size(success_mask, 1) /= size(names, 1))) then
      call log_event( 'Arguments "names" and "success_mask" to function' &
                      // '"ensure_configuration" are different shapes',  &
                      LOG_LEVEL_ERROR )
    end if

    ensure_configuration = .True.

    name_loop: do i = 1, size(names)
      select case(trim( names(i) ))
{ensureSection}
        case default
          write( log_scratch_space, '(A, A, A)' )          &
               "Tried to ensure unrecognised namelist """, &
               trim(names(i)),                             &
               """ was loaded"
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end select

      ensure_configuration = ensure_configuration .and. configuration_found

      if (present(success_mask)) success_mask(i) = configuration_found

    end do name_loop

  end function ensure_configuration

end module configuration_mod
        '''.strip().replace( '~', ' ' )

        firstUsage = '''
~~use foo_config_mod, only : read_foo_namelist, foo_is_loadable, foo_is_loaded
        '''.strip().replace( '~', ' ' )

        firstReadSection = '''
~~~~~~~~~~case ('foo')
            if (foo_is_loadable()) then
              backspace( unit )
              call read_foo_namelist( unit )
            else
              write( log_scratch_space, '(A, A, A)' ) &
                   "Namelist """,                     &
                   trim(buffer),                      &
                   """ can not be read. Too many instances?"
              call log_event( log_scratch_space, LOG_LEVEL_ERROR )
            end if
        '''.strip().replace( '~', ' ' )

        firstEnsureSection = '''
~~~~~~~~case ('foo')
          configuration_found = foo_is_loaded()
        '''.strip().replace( '~', ' ' )

        outputFile = StringIO.StringIO()
        uut = loader.ConfigurationLoader()
        uut.addNamelist( 'foo' )
        uut.writeModule( outputFile )

        inserts = {'usage'         : firstUsage,         \
                   'readSection'   : firstReadSection,   \
                   'ensureSection' : firstEnsureSection}
        self.assertMultiLineEqual( expectedSource.format( **inserts ), \
                                   outputFile.getvalue() )
