#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
# Process previously analysed dependency database. For fun and profit!

from __future__ import print_function;

import os.path

from utilities.path import replaceExtension

###############################################################################
# Process dependency database.
#
class FortranProcessor():
    ###########################################################################
    # Constructor.
    #
    # Arguments:
    #   logger    - A Logger object to write output to.
    #   database  - FortranDatabase object holding details.
    #   objectdir - The directory which holds .o files.
    #   moduledir - The directory which holds .mod files.
    #
    def __init__( self, logger, database, objectDirectory, moduleDirectory):
        self._logger          = logger
        self._database        = database
        self._objectDirectory = objectDirectory
        self._moduleDirectory = moduleDirectory

    ###########################################################################
    # Examine the program unit dependecies and work out the file dependencies.
    #
    # Arguments:
    #   fileStore - FileDependencies object to accept computed dependencies.
    #
    def determineFileDependencies( self, fileStore ):
        currentUnit   = None
        currentObject = None
        for unit, unitFilename, prerequisite, prerequisiteFilename \
                in self._database.getAllDependencies():
            objectFilename = replaceExtension( unitFilename, 'o' )
            objectPathname = os.path.join( self._objectDirectory, \
                                           objectFilename )

            moduleFilename = replaceExtension( prerequisiteFilename, 'mod' )
            modulePathname = os.path.join( self._moduleDirectory, \
                                           moduleFilename )

            fileStore.addFileDependency( objectPathname, modulePathname )

            if unitFilename != currentUnit:
                if currentUnit:
                    fileStore.addFileDependency( currentObject, \
                                                 currentUnit )
                currentUnit = unitFilename
                currentObject = objectPathname

        if currentUnit:
            fileStore.addFileDependency( currentObject, currentUnit )

    ###########################################################################
    # Determine all program units needed to build each program.
    #
    # TODO: Once we have a more recent version of SQLite we could consider
    # doing this at the database level.
    #
    def determineProgramDependencies( self ):
        for program in self._database.getPrograms():
            prerequisites = set()
            self._descend( program, prerequisites )
            unit, unit_file, prereq, prereq_file = self._database.getDependencies( program ).next()
            program_object_file = replaceExtension( unit_file, 'o' )
            prerequisites.add( program_object_file )
            yield program, prerequisites

    ###########################################################################
    def _descend( self, programUnit, prerequisites ):
        for unit, unit_file, prereq, prereq_file \
                              in self._database.getDependencies( programUnit ):
            prereq_object_file = replaceExtension( prereq_file, 'o' )

            self._descend( prereq, prerequisites )
            prerequisites.add( prereq_object_file )
