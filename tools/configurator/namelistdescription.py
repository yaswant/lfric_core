#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
'''
Turns namelist descriptions into namelist modules.
'''

from __future__ import print_function

import pyparsing as parsing

###############################################################################
class NamelistDescriptionException(Exception):
    pass

###############################################################################
class NamelistDescription():
    _moduleTemplate = '''
!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
!> Manages the {listname} namelist.
!>
module {listname}_config_mod

  use constants_mod, only : {kindlist}
  use log_mod,       only : log_event, log_scratch_space, LOG_LEVEL_ERROR

  implicit none

  private
  public :: {publics}

{enumerations}

{variables}

{enumerationHelpers}

  logical :: namelist_loaded = .False.

contains

{enumerationKeyFunctions}

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> \\param [in] file_unit Unit number of the file to read from.
  !>
  subroutine read_{listname}_namelist( file_unit )

{constants}

    implicit none

    integer(i_native), intent(in) :: file_unit

{enumerationKeys}

{namelist}

    integer(i_native) :: condition

    read( file_unit, nml={listname}, iostat=condition, iomsg=log_scratch_space )
    if (condition /= 0) then
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

{interpretEnumerationKey}

{initialisation}

    namelist_loaded = .True.

  end subroutine read_{listname}_namelist

  !> Can this namelist be loaded?
  !>
  !> \\return True if it is possible to load the namelist.
  !>
  function {listname}_is_loadable()

    implicit none

    logical :: {listname}_is_loadable

    {listname}_is_loadable = .not. namelist_loaded

  end function {listname}_is_loadable

  !> Has this namelist been loaded?
  !>
  !> \\return True if the namelist has been loaded.
  !>
  function {listname}_is_loaded()

    implicit none

    logical :: {listname}_is_loaded

    {listname}_is_loaded = namelist_loaded

  end function {listname}_is_loaded

end module {listname}_config_mod
    '''.strip()

    _enumerationTemplate = '  integer(i_native), public, parameter :: {listname}_{enumeration}_{identifier} = {value}'

    _variableTemplate = '  {type}({kind}), public, protected :: {name}'
    _enumVariableTemplate = '  {type}({kind}), public, protected :: {name}'

    _enumerationHelperTemplate = '''  character(str_short), parameter :: {enumeration}_key({key_count}) = (/{keys}/)
  integer(i_native), public, protected :: module_{enumeration}
  equivalence( module_{enumeration}, {enumeration} )'''

    _enumerationKeyTemplate = '    character(str_short) :: {enumeration}'

    _namelistTemplate = '    namelist /{listname}/ {variables}'

    _enumFromKeyTemplate = '''  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> \param[in] key Enumeration key.
  !>
  integer(i_native) function {enumeration}_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    key_index = 1
    do
      if (trim({enumeration}_key(key_index)) == trim(key)) then
        {enumeration}_from_key = key_index + {listname}_{enumeration}_{first_key} - 1
        return
      else
        key_index = key_index + 1
        if (key_index > ubound({enumeration}_key, 1)) then
          write( log_scratch_space, '("Key ''", A, "'' not recognised for {listname} {enumeration}")' ) key
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function {enumeration}_from_key
'''

    _keyFromEnumTemplate = '''  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> \param[in] value Enumeration value.
  !>
  character(str_short) function key_from_{enumeration}( value )

    implicit none

    integer(i_native), intent(in) :: value

    integer(i_native) :: key_index

    key_index = value - {listname}_{enumeration}_{first_key} + 1
    if (key_index < lbound({enumeration}_key, 1) &
        .or. key_index > ubound({enumeration}_key, 1)) then
      write( log_scratch_space, '("Value ", I0, " is not in {listname} {enumeration}")' ) value
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    key_from_{enumeration} = {enumeration}_key( key_index )

  end function key_from_{enumeration}'''

    _interpretEnumerationKeyTemplate = '''    module_{enumeration} = {enumeration}_from_key( {enumeration} )'''

    _constantsTemplate = '    use constants_mod, only : {constants}'

    _initialisationTemplate = '    {name} = {code}'

    _enumerationType = 'integer'
    _enumerationKind = 'native'

    _readableTypeMap = { 'logical' : 'logical', \
                         'integer' : 'integer', \
                         'real'    : 'real',    \
                         'string'  : 'string' }
    _readableKindMap = { 'default' : 'default', \
                         'native'  : 'native',  \
                         'short'   : 'short',   \
                         'long'    : 'long',    \
                         'single'  : 'single',  \
                         'double'  : 'double' }

    _fortranTypeMap = { 'logical' : 'logical', \
                        'integer' : 'integer', \
                        'real'    : 'real',    \
                        'string'  : 'character' }
    _fortranKindMap = { 'default' : 'def',    \
                        'native'  : 'native', \
                        'short'   : 'short',  \
                        'long'    : 'long',   \
                        'single'  : 'single', \
                        'double'  : 'double' }

    _stringLengthMap = { 'default'  : 'str_def', \
                         'filename' : 'str_max_filename' }

    ###########################################################################
    def __init__( self, name ):
        self._name         = name
        self._parameters   = {}
        self._enumerations = {}
        self._computed     = {}
        self._constants    = set()

    ###########################################################################
    def getName( self ):
        return self._name

    ###########################################################################
    def asDict( self ):
        if len(self._parameters) + len(self._enumerations) == 0:
            return {}
        else:
            representation = {}

            for xtype, remainder in self._parameters.items():
                for kind, tail in remainder.items():
                    for name, args in tail.items():
                        representation[name] = [xtype, kind]
                        representation[name].extend( args )

            for name, identifiers in self._enumerations.items():
                representation[name] = ['enumeration', None]
                representation[name].extend( identifiers )

            for name in self._constants:
                representation[name] = ['constant']

            return {self._name : representation}

    ###########################################################################
    def addParameter( self, name, xtype, kind='default', args=[] ):
        if xtype == 'constant':
            self._constants.add( name )
        else:
            if xtype == 'enumeration':
                self._enumerations[name] = args
                xtype = 'integer'
                kind  = 'native'
            elif args:
                self._computed[name] = args

            if not kind:
                kind = 'default'

            if xtype not in self._parameters:
                self._parameters[xtype] = {}
            if kind not in self._parameters[xtype]:
                self._parameters[xtype][kind] = {}
            self._parameters[xtype][kind][name] = args

    ###########################################################################
    def getModuleName ( self ):
        return self._name + '_config_mod'

    ###########################################################################
    def writeModule( self, fileObject ):
        if len(self._parameters) + len(self._enumerations) == 0:
            message = 'Namelist description contains no variables.'
            raise NamelistDescriptionException( message )

        kindset             = {'i_native'}
        variables           = []
        publics             = []
        enumerations        = []
        enumerationHelpers  = []
        enumerationFuncts   = []
        enumerationKeys     = []
        enumInterpreters    = []
        names               = []
        interfaceProcedures = []
        definitions         = {}
        initialisation      = []

        if self._enumerations:
            kindset.add( 'str_short' )

        publics.extend( ['read_{}_namelist'.format( self._name.lower() ), \
                         '{}_is_loadable'.format( self._name.lower() ),   \
                         '{}_is_loaded'.format( self._name.lower() )] )
        evalue = 100
        for xtype, kindBag in self._parameters.items():
            for kind, nameBag in kindBag.items():
                if xtype == 'string':
                    # Strings are a special case because they have a lengh, not
                    # a kind.
                    readableType = kind + '_' + self._readableTypeMap[xtype]
                    fortranType = self._fortranTypeMap[xtype] \
                                  + '(' + self._stringLengthMap[kind] + ')'
                    kindset.add( self._stringLengthMap[kind] )
                else:
                    readableType = self._readableKindMap[kind] \
                                   + '_' + self._readableTypeMap[xtype]
                    fortranKind = self._fortranTypeMap[xtype][0] \
                                  + '_' + self._fortranKindMap[kind]
                    fortranType = self._fortranTypeMap[xtype] \
                                  + '(' + fortranKind + ')'
                    kindset.add( fortranKind )
                humanReadableType = readableType.replace( '_', ' ' )

                for name, args in nameBag.items():
                    names.append( name )

                    # Create the list of variables
                    if xtype == 'string' :
                        inserts = {'type' : self._fortranTypeMap[xtype], \
                                   'kind' : self._stringLengthMap[kind], \
                                   'name' : name}
                    else:
                        inserts = {'type' : self._fortranTypeMap[xtype],      \
                                 'kind' : self._fortranTypeMap[xtype][0]      \
                                          + '_' + self._fortranKindMap[kind], \
                                 'name' : name}
                    if name in self._enumerations:
                        text = self._enumVariableTemplate.format( **inserts )
                    else:
                        text = self._variableTemplate.format( **inserts )
                    variables.append( text )

                    if name in self._enumerations:
                        publics.append( '{}_from_key'.format( name.lower() ) )
                        publics.append( 'key_from_{}'.format( name.lower() ) )

                        inserts = {'enumeration' : name.lower()}
                        text = self._enumerationKeyTemplate.format( **inserts )
                        enumerationKeys.append( text )

                        firstKey = None
                        keys = []
                        for eid in args:
                            if not firstKey: firstKey = eid.lower()
                            keys.append( "'" + eid.lower() + "'" )
                            inserts = {'listname'    : self._name.upper(),  \
                                       'enumeration' : name.upper(),        \
                                       'identifier'  : eid.upper(),         \
                                       'value'       : evalue}
                            text = self._enumerationTemplate.format( **inserts )
                            enumerations.append( text )
                            evalue += 1
                        lastKey = eid.lower()

                        inserts = {'listname'    : self._name.lower(), \
                                   'enumeration' : name.lower()}
                        text = self._interpretEnumerationKeyTemplate.format( **inserts )
                        enumInterpreters.append( text )

                        inserts = {'enumeration' : name.lower(),       \
                                   'listname'    : self._name.lower(), \
                                   'enumeration' : name.lower(),       \
                                   'first_key'   : firstKey,           \
                                   'last_key'    : lastKey,            \
                                   'keys'        : ', '.join( keys ),  \
                                   'key_count'   : len(keys)}
                        text = self._enumerationHelperTemplate.format( **inserts )
                        enumerationHelpers.append( text )

                        text = self._enumFromKeyTemplate.format( **inserts )
                        enumerationFuncts.append( text )
                        text = self._keyFromEnumTemplate.format( **inserts )
                        enumerationFuncts.append( text )
                    elif name in self._computed:
                        if args:
                            inserts = {'name' : name, 'code' : args[0]}
                            code = self._initialisationTemplate.format( **inserts )
                            initialisation.append( code )

        if len(names) > 0:
            noncomputed = list( set(names) - set(self._computed.keys()) )
            noncomputed.sort()
            inserts = {'listname' : self._name, \
                       'variables' : ', '.join(noncomputed)}
            namelist = self._namelistTemplate.format( **inserts )
        else:
            namelist = ''

        if len(self._constants) > 0:
            inserts = {'constants' : ', '.join(self._constants)}
            constants = self._constantsTemplate.format(**inserts)
        else:
            constants = ''

        inserts = {'listname'                : self._name,                    \
                   'kindlist'                : ', '.join( sorted( kindset ) ),\
                   'publics'                 : ', '.join( sorted( publics ) ),\
                   'enumerations'            : '\n'.join(enumerations),       \
                   'variables'               : '\n'.join(sorted(variables)),  \
                   'enumerationHelpers'      : '\n'.join(enumerationHelpers), \
                   'enumerationKeyFunctions' : '\n'.join(enumerationFuncts),  \
                   'enumerationKeys'         : '\n'.join(enumerationKeys),    \
                   'namelist'                : namelist,                      \
                   'interpretEnumerationKey' : '\n'.join(enumInterpreters),   \
                   'constants'               : constants,                     \
                   'initialisation'          : '\n'.join(initialisation)}
        print( self._moduleTemplate.format(**inserts ), end='', \
               file=fileObject )

###############################################################################
class NamelistDescriptionParser():
    '''
    Syntax of namelist description file:

    namelistname ::= alpha[alphanum]*
    file ::= "namelist" namelistname "end" "namelist" namelistname
    '''

    ###########################################################################
    def __init__( self ):
        exclam      = parsing.Literal('!')
        colon       = parsing.Literal(':').suppress()
        openParen   = parsing.Literal('(').suppress()
        closeParen  = parsing.Literal(')').suppress()
        openSquare  = parsing.Literal('[').suppress()
        closeSquare = parsing.Literal(']').suppress()
        comma       = parsing.Literal(',').suppress()

        namelistKeyword = parsing.Literal('namelist').suppress()
        endKeyword      = parsing.Literal('end').suppress()
        typeKeywords \
                  = parsing.oneOf( 'logical integer real string enumeration constant', \
                                   caseless=True )
        kindKeywords \
                   = parsing.oneOf('default native short long single double', \
                                   caseless=True)
        stringKeywords = parsing.oneOf('default, filename', caseless=True)

        name = parsing.Forward()
        def catchName( tokens ):
            name << parsing.oneOf(tokens.asList())
        nameLabel = parsing.Word( parsing.alphas+"_", parsing.alphanums+"_" )
        nameLabel.setParseAction( catchName )
        definitionStart = namelistKeyword + nameLabel
        definitionEnd = endKeyword + namelistKeyword - name.suppress()

        kinddef = kindKeywords ^ stringKeywords
        typedef = typeKeywords('xtype') \
                  + parsing.Optional( openParen \
                                      - kinddef('xkind') \
                                      - closeParen )

        label = parsing.Word( parsing.alphas+"_", parsing.alphanums+"_" )
        argument = label ^ parsing.quotedString.addParseAction(parsing.removeQuotes)
        argumentList = parsing.Group( openSquare + argument \
                                      + parsing.ZeroOrMore( comma + argument)
                                      + closeSquare ).setResultsName('xargs')

        variable = parsing.Group( label + colon + typedef \
                                  + parsing.Optional( argumentList ) )

        definition = parsing.Group( definitionStart                           \
                              + parsing.Dict( parsing.OneOrMore( variable ) ) \
                              + definitionEnd )

        self._parser = parsing.Dict( definition )

        comment = exclam - parsing.restOfLine
        self._parser.ignore( comment )

    ###########################################################################
    def parseFile ( self, fileObject ):
        try:
            parseTree = self._parser.parseFile( fileObject )
        except (parsing.ParseException, parsing.ParseSyntaxException) as err:
            message = '\n{}\n{}\n{}'.format( err.line, \
                                           " "*(err.column-1) + "^", \
                                           err )
            raise NamelistDescriptionException( message )

        result = []
        for name, variables in parseTree.items():
            description = NamelistDescription( name )

            for key, value in variables.items():
                if isinstance(value, parsing.ParseResults) :
                    description.addParameter( key,         \
                                              value.xtype, \
                                              value.xkind, \
                                              value.xargs )
                else:
                    description.addParameter( key, value )

            result.append( description )

        return result