#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
'''
Implements a Jinja2 filter which removes duplicates from a Cylc task schedule.
'''
import re

def deduplicateSchedule( schedule ):
    '''
    Takes a Cylc task schedule and removes duplicates.

    @param [in] schedule String containing Cylc task schedule.
    @return Deduplicated version of schedule.
    '''
    dependencyPattern = re.compile( r'\s*(\S+)\s*=>\s*(\S+)' )
    dependencySet = set()
    for match in dependencyPattern.finditer( schedule ):
        dependencySet.add( (match.group(1), match.group(2)) )

    newSchedule = ['{indent}{prerequisite} => {result}'.format(indent=' '*12, \
                                                          prerequisite=first, \
                                                          result=second) \
                   for first, second in dependencySet]
    return '\n'.join( newSchedule )
