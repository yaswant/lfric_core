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
    Takes a Cylc task schedule and removes duplicates.  It also combines
    scheduling dependencies into groups with the same cycling pattern.

    For example, if schedule is:
       a => b
       e => f
       [[[R1]]]
           graph = """
           a => b
           c => d
           """
       [[[R3/P1/3]]]
           graph = """
           a[-P1] => a
           """
     then the returned schedule would be:
        [[[R1]]]
           graph = """
              a => b
              c => d
              e => f
              """
        [[[R3/P1/3]]]
           graph = """
              a[-P1] => a
              """
    i.e., tasks at the beginning of the schedule that aren't explicitly associated
    with a cycle period are put into the [[[R1]]] group of tasks. 

    @param [in] schedule String containing Cylc task schedule.
    @return Deduplicated version of schedule.
    '''
    cyclePattern = re.compile( r'\s*(\[\[\[\S*\]\]\])(.*?)(?=\[\[\[)', re.S )
    dependencyPattern = re.compile( r'\s*(\S+)\s*=>\s*(\S+)' )

    scheduleDict={}
    for match in cyclePattern.finditer( '[[[]]]' + schedule + '[[['):
        k,v = match.group(1), match.group(2)
        if k == '[[[]]]': k='[[[R1]]]'
        if k in scheduleDict.keys():
            scheduleDict[k] = scheduleDict[k] + v
        else:
            scheduleDict[k] =  v

    newSchedule=''
    for k,v in scheduleDict.items():
        dependencySet = set()
        for match in dependencyPattern.finditer( v ):
            dependencySet.add( (match.group(1), match.group(2)) )

        scheduleDict[k] = ['{indent}{prerequisite} => {result}'.format(indent=' '*12, \
                                                          prerequisite=first, \
                                                          result=second) \
                               for first, second in dependencySet]
        newSchedule = '{old}\n{indent1}{cyclegroup}\n{indent2}graph="""\n{schedulegraph}\n{indent2}"""\n'.format( \
                old=newSchedule, \
                indent1=' '*8,  \
                indent2=' '*12, \
                cyclegroup=k, \
                schedulegraph='\n'.join(scheduleDict[k]) \
                )

    return newSchedule



