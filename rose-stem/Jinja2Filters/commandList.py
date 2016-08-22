#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
'''
Implements a Jinja2 filter which takes a number of strings and joins them with
the filter text.

e.g.

' foo ' | commandList( 'comm1', ['comm2, 'comm3'] )
=> comm1 foo comm2 foo comm3
'''
def commandList( separator, *fragments ):
    '''
    Takes an arbitrary collection of strings and joins them with the filter
    text into a single string.

    @param [in] separator Filter text string
    @param [in] fragments Multiple strings or lists of strings.
    '''
    commands = []
    for fragment in fragments:
        if hasattr(fragment, '__iter__'):
            commands.extend( fragment )
        else:
            commands.append( fragment )
    return separator.join( commands )
