#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
###############################################################################
# A library of logging tools.

from __future__ import print_function

import abc

###############################################################################
# Abstract logging provider class
class Logger():
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def logEvent( self, description ):
        raise NotImplementedError( 'Logger.logEvent not implemented' )

###############################################################################
# Silent (null) logger
class NoLogger( Logger ):
    def logEvent( self, description ):
        pass

###############################################################################
# Print to file logger
class PrintLogger( Logger ):
    def __init__(self, stream ):
        self.stream = stream

    def logEvent( self, description ):
        print( '{}'.format( description ), file=self.stream )
