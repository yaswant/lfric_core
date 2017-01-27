#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
'''
Implements a Jinja2 filter to run a macro specified by a string.
'''
from jinja2 import contextfilter
import re

@contextfilter
def executeMacroCrun(context, call):
    '''
    Takes a string and executes it as though it were a Jinja2 macro call, but 
    overrides keyword do_crun to be True

    The call string has the syntax <macro name>([<argument>]...).

    Arguments can be either position or keyword based.

    @param [inout] context Jinja2 instance to run macro against.
    @param [in]    call    Invokation string.
    @return String resulting from calling the macro.
    '''
    if call.find('(') == -1:
        macroName = call
        arguments = ''
    else:
        macroName = call[:call.index('(')]
        arguments = re.split(', *', call[call.index('(')+1:call.rindex(')')])

    normalArguments  = [argument for argument in arguments \
                        if argument.find('=') == -1]
    keywordArguments = [argument for argument in arguments \
                        if argument.find('=') != -1]

    argumentList = []
    for argument in normalArguments:
        if argument[0] == '"':
            argumentList.append( argument[1:-2] )
        else:
            argumentList.append( argument )

    argumentDictionary = {}
    for argument in keywordArguments:
        key, value = re.split(' *= *', argument)
        argumentDictionary[key] = value

    argumentDictionary['do_crun'] = True
    return context.vars[macroName]( *argumentList, **argumentDictionary )
