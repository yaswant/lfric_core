#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

from __future__ import print_function

import jinja2    as jinja

import jinjamacros

##############################################################################
class NamelistFeigner():
    def __init__( self ):
        self._engine = jinja.Environment( \
                   loader=jinja.PackageLoader( 'configurator', 'templates'), \
                   extensions=['jinja2.ext.do'] )
        self._engine.filters['decorate'] = jinjamacros.decorateMacro

        self._namelists = {}
        self._kinds     = set( ['i_native'] )
        self._arguments = {}

    def addNamelist( self, namelists ):
        for item in namelists:
            name = item.getName()
            self._namelists[name] = item

            self._arguments[name] = []
            for param, fortranType in item.getParameters().iteritems():
                if param not in item.getComputations():
                    self._arguments[name].append( param )
                self._kinds.add( fortranType.kind )

    def writeModule( self, moduleFile ):
        inserts = { 'namelists' : self._namelists, \
                    'kinds'     : self._kinds,     \
                    'arguments' : self._arguments }

        template = self._engine.get_template( 'feign_config.f90' )
        print( template.render( inserts ), file=moduleFile )
