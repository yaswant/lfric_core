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

import configurator.namelistdescription as namelist
import configurator.namelistfeigner     as feigner

###############################################################################
class FeignerTest( unittest.TestCase ):
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
! Handles the feigning of namelists.
!
module feign_config_mod

  use constants_mod, only : i_native

  implicit none

  private
  public :: 

  integer(i_native), parameter :: temporary_unit = 3

contains

end module feign_config_mod
        '''.strip().replace( '~', ' ' )

        outputFile = StringIO.StringIO()
        uut = feigner.NamelistFeigner()
        uut.writeModule( outputFile )

        self.assertMultiLineEqual( expectedSource + '\n', \
                                   outputFile.getvalue() )

    ###########################################################################
    def testSimple( self ):
        expectedSource = '''
!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
! Handles the feigning of namelists.
!
module feign_config_mod

  use constants_mod, only : i_def, i_native, l_def, r_double, str_def

  implicit none

  private
  public :: feign_simple_config

  integer(i_native), parameter :: temporary_unit = 3

contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine feign_simple_config( foo, &
                                  bar, &
                                  baz, &
                                  fred )

    use simple_config_mod, only : read_simple_namelist

    implicit none

    integer(i_def), intent(in) :: foo
    real(r_double), intent(in) :: bar
    character(*), intent(in) :: baz
    logical(l_def), intent(in) :: fred

    integer(i_native) :: condition

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) then
      write( 6, '("feign_simple_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&simple")' )
    write( temporary_unit, '("foo = ", I0)' ) foo
    write( temporary_unit, '("bar = ", E14.7)' ) bar
    write( temporary_unit, '("baz = ''", A, "''")' ) baz
    write( temporary_unit, '("fred = ", L)' ) fred
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_simple_namelist( temporary_unit )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop 'feign_simple_config: Unable to close temporary file'

  end subroutine feign_simple_config

end module feign_config_mod
        '''.strip().replace( '~', ' ' )

        simple = namelist.NamelistDescription( 'simple' )
        simple.addParameter( 'foo', 'integer', 'default' )
        simple.addParameter( 'bar', 'real', 'double' )
        simple.addParameter( 'baz', 'string' )
        simple.addParameter( 'qux', 'constant' )
        simple.addParameter( 'fred', 'logical' )

        outputFile = StringIO.StringIO()
        uut = feigner.NamelistFeigner()
        uut.addNamelist( [simple] )
        uut.writeModule( outputFile )

        self.assertMultiLineEqual( expectedSource + '\n', \
                                   outputFile.getvalue() )

    ###########################################################################
    def testEnumeration( self ):
        expectedSource = '''
!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
! Handles the feigning of namelists.
!
module feign_config_mod

  use constants_mod, only : i_native

  implicit none

  private
  public :: feign_enum_config

  integer(i_native), parameter :: temporary_unit = 3

contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine feign_enum_config( thing )

    use enum_config_mod, only : read_enum_namelist, &
                                key_from_thing, &
                                thing_from_key

    implicit none

    integer(i_native), intent(in) :: thing

    integer(i_native) :: condition

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) then
      write( 6, '("feign_enum_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&enum")' )
    write( temporary_unit, '("thing = ''", A, "''")' ) key_from_thing( thing )
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_enum_namelist( temporary_unit )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop 'feign_enum_config: Unable to close temporary file'

  end subroutine feign_enum_config

end module feign_config_mod
        '''.strip().replace( '~', ' ' )

        enumable = namelist.NamelistDescription( 'enum' )
        enumable.addParameter( 'thing', 'enumeration', args=['one', 'two'] )

        outputFile = StringIO.StringIO()
        uut = feigner.NamelistFeigner()
        uut.addNamelist( [enumable] )
        uut.writeModule( outputFile )

        self.assertMultiLineEqual( expectedSource + '\n', \
                                   outputFile.getvalue() )

    ###########################################################################
    def testComputed( self ):
        expectedSource = '''
!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
! Handles the feigning of namelists.
!
module feign_config_mod

  use constants_mod, only : i_def, i_native

  implicit none

  private
  public :: feign_computed_config

  integer(i_native), parameter :: temporary_unit = 3

contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine feign_computed_config( teapot, &
                                    cheese )

    use computed_config_mod, only : read_computed_namelist

    implicit none

    integer(i_def), intent(in) :: teapot
    integer(i_def), intent(in) :: cheese

    integer(i_native) :: condition

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) then
      write( 6, '("feign_computed_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&computed")' )
    write( temporary_unit, '("teapot = ", I0)' ) teapot
    write( temporary_unit, '("cheese = ", I0)' ) cheese
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_computed_namelist( temporary_unit )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop 'feign_computed_config: Unable to close temporary file'

  end subroutine feign_computed_config

end module feign_config_mod
        '''.strip().replace( '~', ' ' )

        simple = namelist.NamelistDescription( 'computed' )
        simple.addParameter( 'teapot', 'integer', 'default' )
        simple.addParameter( 'cheese', 'integer', 'default' )
        simple.addParameter( 'biscuits', 'integer', 'default', ['teapot + cheese'] )

        outputFile = StringIO.StringIO()
        uut = feigner.NamelistFeigner()
        uut.addNamelist( [simple] )
        uut.writeModule( outputFile )

        self.assertMultiLineEqual( expectedSource + '\n', \
                                   outputFile.getvalue() )

    ###########################################################################
    def testEverything( self ):
        expectedSource = '''
!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
! Handles the feigning of namelists.
!
module feign_config_mod

  use constants_mod, only : i_native, l_def, r_def, str_max_filename

  implicit none

  private
  public :: feign_everything_config

  integer(i_native), parameter :: temporary_unit = 3

contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine feign_everything_config( cake, &
                                      teapot, &
                                      cheese, &
                                      fish, &
                                      tail )

    use everything_config_mod, only : read_everything_namelist, &
                                      key_from_teapot, &
                                      teapot_from_key

    implicit none

    character(*), intent(in) :: cake
    integer(i_native), intent(in) :: teapot
    logical(l_def), intent(in) :: cheese
    real(r_def), intent(in) :: fish
    integer(i_native), intent(in) :: tail

    integer(i_native) :: condition

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) then
      write( 6, '("feign_everything_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&everything")' )
    write( temporary_unit, '("cake = ''", A, "''")' ) cake
    write( temporary_unit, '("teapot = ''", A, "''")' ) key_from_teapot( teapot )
    write( temporary_unit, '("cheese = ", L)' ) cheese
    write( temporary_unit, '("fish = ", E14.7)' ) fish
    write( temporary_unit, '("tail = ", I0)' ) tail
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_everything_namelist( temporary_unit )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop 'feign_everything_config: Unable to close temporary file'

  end subroutine feign_everything_config

end module feign_config_mod
        '''.strip().replace( '~', ' ' )

        everything = namelist.NamelistDescription( 'everything' )
        everything.addParameter( 'cake', 'string', 'filename' )
        everything.addParameter( 'teapot', 'enumeration', args=['brown', \
                                                                'steel'] )
        everything.addParameter( 'cheese', 'logical' )
        everything.addParameter( 'fish', 'real' )
        everything.addParameter( 'wibble', 'constant' )
        everything.addParameter( 'yarn', 'real', 'default', \
                                 ['fish * wibble / 180.0_r_def'] )
        everything.addParameter( 'tail', 'integer', 'native' )

        outputFile = StringIO.StringIO()
        uut = feigner.NamelistFeigner()
        uut.addNamelist( [everything] )
        uut.writeModule( outputFile )

        self.assertMultiLineEqual( expectedSource + '\n', \
                                   outputFile.getvalue() )
