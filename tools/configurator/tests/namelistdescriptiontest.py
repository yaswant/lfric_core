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

import configurator.namelistdescription as description

###############################################################################
class NamelistDescriptionTest( unittest.TestCase ):
    ###########################################################################
    def setUp( self ):
        self.maxDiff = None

    ###########################################################################
    def tearDown( self ):
        pass

    ###########################################################################
    def testModuleWriteEmpty( self ):
        outputFile = StringIO.StringIO()

        uut = description.NamelistDescription( 'test' )
        self.assertRaises( description.NamelistDescriptionException, \
                           uut.writeModule, outputFile )

    ###########################################################################
    def testModuleWriteOneOfEach( self ):
        expectedSource = '''
!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
!> Manages the test namelist.
!>
module test_config_mod

  use constants_mod, only : i_def, i_long, i_native, i_short, r_def, r_double, r_single, str_def, str_max_filename, str_short
  use log_mod,       only : log_event, log_scratch_space, LOG_LEVEL_ERROR

  implicit none

  private
  public :: enum_from_key, key_from_enum, read_test_namelist, test_is_loadable, test_is_loaded

  integer(i_native), public, parameter :: TEST_ENUM_ONE = 100
  integer(i_native), public, parameter :: TEST_ENUM_TWO = 101
  integer(i_native), public, parameter :: TEST_ENUM_THREE = 102

  character(str_def), public, protected :: dstr
  character(str_def), public, protected :: vstr
  character(str_max_filename), public, protected :: fstr
  integer(i_def), public, protected :: dint
  integer(i_def), public, protected :: vint
  integer(i_long), public, protected :: lint
  integer(i_native), public, protected :: enum
  integer(i_short), public, protected :: sint
  real(r_def), public, protected :: dreal
  real(r_def), public, protected :: vreal
  real(r_double), public, protected :: lreal
  real(r_single), public, protected :: sreal

  character(str_short), parameter :: enum_key(3) = (/'one', 'two', 'three'/)
  integer(i_native), public, protected :: module_enum
  equivalence( module_enum, enum )

  logical :: namelist_loaded = .False.

contains

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> \param[in] key Enumeration key.
  !>
  integer(i_native) function enum_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    key_index = 1
    do
      if (trim(enum_key(key_index)) == trim(key)) then
        enum_from_key = key_index + test_enum_one - 1
        return
      else
        key_index = key_index + 1
        if (key_index > ubound(enum_key, 1)) then
          write( log_scratch_space, '("Key ''", A, "'' not recognised for test enum")' ) key
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function enum_from_key

  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> \param[in] value Enumeration value.
  !>
  character(str_short) function key_from_enum( value )

    implicit none

    integer(i_native), intent(in) :: value

    integer(i_native) :: key_index

    key_index = value - test_enum_one + 1
    if (key_index < lbound(enum_key, 1) &
        .or. key_index > ubound(enum_key, 1)) then
      write( log_scratch_space, '("Value ", I0, " is not in test enum")' ) value
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    key_from_enum = enum_key( key_index )

  end function key_from_enum

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> \param [in] file_unit Unit number of the file to read from.
  !>
  subroutine read_test_namelist( file_unit )



    implicit none

    integer(i_native), intent(in) :: file_unit

    character(str_short) :: enum

    namelist /test/ dint, dreal, dstr, enum, fstr, lint, lreal, sint, sreal, vint, vreal, vstr

    integer(i_native) :: condition

    read( file_unit, nml=test, iostat=condition, iomsg=log_scratch_space )
    if (condition /= 0) then
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    module_enum = enum_from_key( enum )



    namelist_loaded = .True.

  end subroutine read_test_namelist

  !> Can this namelist be loaded?
  !>
  !> \\return True if it is possible to load the namelist.
  !>
  function test_is_loadable()

    implicit none

    logical :: test_is_loadable

    test_is_loadable = .not. namelist_loaded

  end function test_is_loadable

  !> Has this namelist been loaded?
  !>
  !> \\return True if the namelist has been loaded.
  !>
  function test_is_loaded()

    implicit none

    logical :: test_is_loaded

    test_is_loaded = namelist_loaded

  end function test_is_loaded

end module test_config_mod
        '''.strip()

        outputFile = StringIO.StringIO()

        uut = description.NamelistDescription( 'test' )
        uut.addParameter( 'vint', 'integer'             )
        uut.addParameter( 'dint', 'integer', 'default'  )
        uut.addParameter( 'sint', 'integer', 'short'    )
        uut.addParameter( 'lint', 'integer', 'long'     )
        uut.addParameter( 'vreal', 'real'               )
        uut.addParameter( 'dreal', 'real',   'default'  )
        uut.addParameter( 'sreal', 'real',   'single'   )
        uut.addParameter( 'lreal', 'real',   'double'   )
        uut.addParameter( 'vstr', 'string'              )
        uut.addParameter( 'dstr', 'string',  'default'  )
        uut.addParameter( 'fstr', 'string',  'filename' )
        uut.addParameter( 'enum', 'enumeration', None, ['one', 'two', 'three'] )
        uut.writeModule( outputFile )

        self.assertMultiLineEqual( expectedSource, outputFile.getvalue() )

    ###########################################################################
    def testModuleWriteGrowing( self ):
        firstExpectedSource = '''
!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
!> Manages the test namelist.
!>
module test_config_mod

  use constants_mod, only : i_def, i_native
  use log_mod,       only : log_event, log_scratch_space, LOG_LEVEL_ERROR

  implicit none

  private
  public :: read_test_namelist, test_is_loadable, test_is_loaded



  integer(i_def), public, protected :: foo



  logical :: namelist_loaded = .False.

contains



  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> \param [in] file_unit Unit number of the file to read from.
  !>
  subroutine read_test_namelist( file_unit )



    implicit none

    integer(i_native), intent(in) :: file_unit



    namelist /test/ foo

    integer(i_native) :: condition

    read( file_unit, nml=test, iostat=condition, iomsg=log_scratch_space )
    if (condition /= 0) then
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if





    namelist_loaded = .True.

  end subroutine read_test_namelist

  !> Can this namelist be loaded?
  !>
  !> \\return True if it is possible to load the namelist.
  !>
  function test_is_loadable()

    implicit none

    logical :: test_is_loadable

    test_is_loadable = .not. namelist_loaded

  end function test_is_loadable

  !> Has this namelist been loaded?
  !>
  !> \\return True if the namelist has been loaded.
  !>
  function test_is_loaded()

    implicit none

    logical :: test_is_loaded

    test_is_loaded = namelist_loaded

  end function test_is_loaded

end module test_config_mod
        '''.strip()

        secondExpectedSource = '''
!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
!> Manages the test namelist.
!>
module test_config_mod

  use constants_mod, only : i_def, i_native, r_def
  use log_mod,       only : log_event, log_scratch_space, LOG_LEVEL_ERROR

  implicit none

  private
  public :: read_test_namelist, test_is_loadable, test_is_loaded



  integer(i_def), public, protected :: foo
  real(r_def), public, protected :: bar



  logical :: namelist_loaded = .False.

contains



  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> \param [in] file_unit Unit number of the file to read from.
  !>
  subroutine read_test_namelist( file_unit )



    implicit none

    integer(i_native), intent(in) :: file_unit



    namelist /test/ bar, foo

    integer(i_native) :: condition

    read( file_unit, nml=test, iostat=condition, iomsg=log_scratch_space )
    if (condition /= 0) then
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if





    namelist_loaded = .True.

  end subroutine read_test_namelist

  !> Can this namelist be loaded?
  !>
  !> \\return True if it is possible to load the namelist.
  !>
  function test_is_loadable()

    implicit none

    logical :: test_is_loadable

    test_is_loadable = .not. namelist_loaded

  end function test_is_loadable

  !> Has this namelist been loaded?
  !>
  !> \\return True if the namelist has been loaded.
  !>
  function test_is_loaded()

    implicit none

    logical :: test_is_loaded

    test_is_loaded = namelist_loaded

  end function test_is_loaded

end module test_config_mod
        '''.strip()

        outputFile = StringIO.StringIO()

        uut = description.NamelistDescription( 'test' )
        uut.addParameter( 'foo', 'integer' )
        uut.writeModule( outputFile )

        self.assertMultiLineEqual( firstExpectedSource, outputFile.getvalue() )

        outputFile = StringIO.StringIO()
        uut.addParameter( 'bar', 'real', 'default' )
        uut.writeModule( outputFile )

        self.assertMultiLineEqual( secondExpectedSource, outputFile.getvalue() )

    def testEnumerationOnly( self ):
        expectedSource = '''
!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
!> Manages the enum namelist.
!>
module enum_config_mod

  use constants_mod, only : i_native, str_short
  use log_mod,       only : log_event, log_scratch_space, LOG_LEVEL_ERROR

  implicit none

  private
  public :: enum_is_loadable, enum_is_loaded, key_from_value, read_enum_namelist, value_from_key

  integer(i_native), public, parameter :: ENUM_VALUE_ONE = 100
  integer(i_native), public, parameter :: ENUM_VALUE_TWO = 101
  integer(i_native), public, parameter :: ENUM_VALUE_THREE = 102

  integer(i_native), public, protected :: value

  character(str_short), parameter :: value_key(3) = (/'one', 'two', 'three'/)
  integer(i_native), public, protected :: module_value
  equivalence( module_value, value )

  logical :: namelist_loaded = .False.

contains

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> \param[in] key Enumeration key.
  !>
  integer(i_native) function value_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    key_index = 1
    do
      if (trim(value_key(key_index)) == trim(key)) then
        value_from_key = key_index + enum_value_one - 1
        return
      else
        key_index = key_index + 1
        if (key_index > ubound(value_key, 1)) then
          write( log_scratch_space, '("Key ''", A, "'' not recognised for enum value")' ) key
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function value_from_key

  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> \param[in] value Enumeration value.
  !>
  character(str_short) function key_from_value( value )

    implicit none

    integer(i_native), intent(in) :: value

    integer(i_native) :: key_index

    key_index = value - enum_value_one + 1
    if (key_index < lbound(value_key, 1) &
        .or. key_index > ubound(value_key, 1)) then
      write( log_scratch_space, '("Value ", I0, " is not in enum value")' ) value
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    key_from_value = value_key( key_index )

  end function key_from_value

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> \param [in] file_unit Unit number of the file to read from.
  !>
  subroutine read_enum_namelist( file_unit )



    implicit none

    integer(i_native), intent(in) :: file_unit

    character(str_short) :: value

    namelist /enum/ value

    integer(i_native) :: condition

    read( file_unit, nml=enum, iostat=condition, iomsg=log_scratch_space )
    if (condition /= 0) then
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    module_value = value_from_key( value )



    namelist_loaded = .True.

  end subroutine read_enum_namelist

  !> Can this namelist be loaded?
  !>
  !> \\return True if it is possible to load the namelist.
  !>
  function enum_is_loadable()

    implicit none

    logical :: enum_is_loadable

    enum_is_loadable = .not. namelist_loaded

  end function enum_is_loadable

  !> Has this namelist been loaded?
  !>
  !> \\return True if the namelist has been loaded.
  !>
  function enum_is_loaded()

    implicit none

    logical :: enum_is_loaded

    enum_is_loaded = namelist_loaded

  end function enum_is_loaded

end module enum_config_mod
        '''.strip()

        outputFile = StringIO.StringIO()

        uut = description.NamelistDescription( 'enum' )
        uut.addParameter( 'value', 'enumeration', None, ['one', 'two', 'three'] )
        uut.writeModule( outputFile )

        self.assertMultiLineEqual( expectedSource, outputFile.getvalue() )

    ###########################################################################
    def testModuleWriteComputed( self ):
        expectedSource = '''
!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
!> Manages the teapot namelist.
!>
module teapot_config_mod

  use constants_mod, only : i_native, r_def
  use log_mod,       only : log_event, log_scratch_space, LOG_LEVEL_ERROR

  implicit none

  private
  public :: read_teapot_namelist, teapot_is_loadable, teapot_is_loaded



  real(r_def), public, protected :: bar
  real(r_def), public, protected :: foo



  logical :: namelist_loaded = .False.

contains



  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> \param [in] file_unit Unit number of the file to read from.
  !>
  subroutine read_teapot_namelist( file_unit )



    implicit none

    integer(i_native), intent(in) :: file_unit



    namelist /teapot/ foo

    integer(i_native) :: condition

    read( file_unit, nml=teapot, iostat=condition, iomsg=log_scratch_space )
    if (condition /= 0) then
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if



    bar = foo ** 2

    namelist_loaded = .True.

  end subroutine read_teapot_namelist

  !> Can this namelist be loaded?
  !>
  !> \\return True if it is possible to load the namelist.
  !>
  function teapot_is_loadable()

    implicit none

    logical :: teapot_is_loadable

    teapot_is_loadable = .not. namelist_loaded

  end function teapot_is_loadable

  !> Has this namelist been loaded?
  !>
  !> \\return True if the namelist has been loaded.
  !>
  function teapot_is_loaded()

    implicit none

    logical :: teapot_is_loaded

    teapot_is_loaded = namelist_loaded

  end function teapot_is_loaded

end module teapot_config_mod
        '''.strip()

        outputFile = StringIO.StringIO()

        uut = description.NamelistDescription( 'teapot' )
        uut.addParameter( 'foo', 'real', 'default' )
        uut.addParameter( 'bar', 'real', 'default', ['foo ** 2'] )
        uut.writeModule( outputFile )

        self.assertMultiLineEqual( expectedSource, outputFile.getvalue() )

    ###########################################################################
    def testModuleWriteConstant( self ):
        expectedSource = '''
!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
!> Manages the cheese namelist.
!>
module cheese_config_mod

  use constants_mod, only : i_native, r_def
  use log_mod,       only : log_event, log_scratch_space, LOG_LEVEL_ERROR

  implicit none

  private
  public :: cheese_is_loadable, cheese_is_loaded, read_cheese_namelist



  real(r_def), public, protected :: fred
  real(r_def), public, protected :: wilma



  logical :: namelist_loaded = .False.

contains



  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> \param [in] file_unit Unit number of the file to read from.
  !>
  subroutine read_cheese_namelist( file_unit )

    use constants_mod, only : FUDGE

    implicit none

    integer(i_native), intent(in) :: file_unit



    namelist /cheese/ fred

    integer(i_native) :: condition

    read( file_unit, nml=cheese, iostat=condition, iomsg=log_scratch_space )
    if (condition /= 0) then
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if



    wilma = fred * FUDGE

    namelist_loaded = .True.

  end subroutine read_cheese_namelist

  !> Can this namelist be loaded?
  !>
  !> \\return True if it is possible to load the namelist.
  !>
  function cheese_is_loadable()

    implicit none

    logical :: cheese_is_loadable

    cheese_is_loadable = .not. namelist_loaded

  end function cheese_is_loadable

  !> Has this namelist been loaded?
  !>
  !> \\return True if the namelist has been loaded.
  !>
  function cheese_is_loaded()

    implicit none

    logical :: cheese_is_loaded

    cheese_is_loaded = namelist_loaded

  end function cheese_is_loaded

end module cheese_config_mod
        '''.strip()

        outputFile = StringIO.StringIO()

        uut = description.NamelistDescription( 'cheese' )
        uut.addParameter( 'FUDGE', 'constant' )
        uut.addParameter( 'fred', 'real', 'default' )
        uut.addParameter( 'wilma', 'real', 'default', ['fred * FUDGE'] )
        uut.writeModule( outputFile )

        self.assertMultiLineEqual( expectedSource, outputFile.getvalue() )

###############################################################################
class NamelistDescriptionParserTest( unittest.TestCase ):
    ###########################################################################
    def setUp( self ):
        pass

    ###########################################################################
    def tearDown( self ):
        pass

    ###########################################################################
    def testParserGoodFile( self ):
        inputFile = StringIO.StringIO('''
        ! Initial comment followed by blank line

        namelist fred
            first_thing : string ! Trailing comment
            second      : integer
            filename    : string(filename)
            choices     : enumeration[ foo, bar, baz, qux ]
        end namelist fred
        ''')

        uut = description.NamelistDescriptionParser()
        result = {}
        for namelist in uut.parseFile( inputFile ):
            result.update( namelist.asDict() )

        self.assertEqual( {'fred'  : {'first_thing' : ['string', 'default'],\
                                      'second' : ['integer', 'default'],    \
                                      'filename' : ['string', 'filename'],  \
                                      'choices' : ['enumeration', None, \
                                               'foo', 'bar', 'baz', 'qux']}}, \
                          result )

    ###########################################################################
    def testOnlyEnumeration( self ):
        inputFile = StringIO.StringIO('''
        namelist barney
            stuff : enumeration[ one, two, three ]
        end namelist barney
        ''')

        uut = description.NamelistDescriptionParser()
        result = {}
        for namelist in uut.parseFile( inputFile ):
            result.update( namelist.asDict() )

        self.assertEqual( {'barney' : {'stuff' : ['enumeration', None, \
                                                  'one', 'two', 'three']}}, \
                          result )

    ###########################################################################
    def testMismatchedNames( self ):
        inputFile = StringIO.StringIO('''
        namelist wilma
            junk : integer
        end namelist betty
        ''')

        uut = description.NamelistDescriptionParser()
        self.assertRaises( description.NamelistDescriptionException, \
                           uut.parseFile, inputFile )

    ###########################################################################
    def testComputedFields( self ):
        inputFile = StringIO.StringIO('''
        namelist teapot
           foo : real
           bar : real(default)['foo ** 2']
           baz : real['PI * foo']
        end namelist teapot
        ''')

        result = {}
        uut = description.NamelistDescriptionParser()
        for namelist in uut.parseFile( inputFile ):
            result.update( namelist.asDict() )

        self.assertEqual({'teapot' : {'foo' : ['real', 'default'], \
                                      'bar' : ['real', 'default', 'foo ** 2'], \
                                      'baz' : ['real', 'default', 'PI * foo']}}, \
                         result)

    ###########################################################################
    def testComputedFields( self ):
        inputFile = StringIO.StringIO('''
        namelist cheese
           FUDGE : constant
           fred  : real
           wilma : real['fred * FUDGE']
        end namelist cheese
        ''')

        result = {}
        uut = description.NamelistDescriptionParser()
        for namelist in uut.parseFile( inputFile ):
            result.update( namelist.asDict() )

        self.assertEqual({'cheese' : {'fred' : ['real', 'default'], \
                                      'FUDGE' : ['constant'], \
                                      'wilma' : ['real', 'default', 'fred * FUDGE']}}, \
                         result)

###############################################################################
if __name__ == '__main__':
    unittest.main()
