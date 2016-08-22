#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
'''
Takes the results of parsing Cylc log files with parser.py and render them in
some attractive fashion.
'''
from __future__ import print_function

from abc import ABCMeta, abstractmethod
import datetime
from jinja2 import Environment, PackageLoader

##############################################################################
# Index rendering
##############################################################################
class IndexRenderer(object):
    '''
    Index page renderers inherit from this class.
    '''
    __MetaClass__ = ABCMeta

    @abstractmethod
    def render( self, stream ):
        '''
        Render to the specified file object.
        '''
        pass

##############################################################################
class HtmlIndexRenderer(IndexRenderer):
    '''
    Render the index page to an HTML document.
    '''
    def __init__( self, indexer, suiteUrl ):
        self._indexer   = indexer
        self._suiteUrl = suiteUrl

    def render( self, stream ):
        variables = {'cronout'       : self._indexer.cronOut,       \
                     'crontimestamp' : self._indexer.cronTimestamp, \
                     'pages'         : self._indexer.pages,         \
                     'subsites'      : self._indexer.subsites,      \
                     'oldthreshold'  : datetime.datetime.utcnow()   \
                                       - datetime.timedelta( hours=48 ), \
                     'suiteurl'      : self._suiteUrl}

        environment = Environment( loader=PackageLoader( 'parserender', \
                                                         'templates' ) )
        template = environment.get_template( 'index.html' )
        print( template.render( variables ), file=stream )

##############################################################################
# Compile output rendering
##############################################################################
class CompileRenderer(object):
    '''
    Compile event renderes inherit from this class.
    '''
    __MetaClass__ = ABCMeta

    @abstractmethod
    def render( self, ignoreCodes, stream ):
        '''
        Render to the specified file object.
        '''
        pass

##############################################################################
class HtmlCompileRenderer(CompileRenderer):
    '''
    Render the results of a compile to an HTML document.
    '''
    def __init__( self, context, outParser, errParser ):
        self._context   = context
        self._outParser = outParser
        self._errParser = errParser

    def render( self, ignoreCodes, stream ):
        if ignoreCodes:
            eventBuckets = self._errParser.getEvents( ignoreCodes=ignoreCodes, \
                                                      groupBy='filename' )
        else:
            eventBuckets = self._errParser.getEvents( groupBy='filename' )

        variables = {'context'      : self._context, \
                     'compiler'     : self._outParser.compiler,  \
                     'timestamp'    : self._outParser.started, \
                     'eventBuckets' : eventBuckets}

        environment = Environment(loader=PackageLoader('parserender', \
                                                       'templates'))
        template = environment.get_template( 'compilelog.html' )
        print( template.render( variables ), file = stream )

##############################################################################
# Run output rendering
##############################################################################
class RunRenderer(object):
    '''
    Run-time event renderers inherit from this class.
    '''
    __MetaClass__ = ABCMeta

    @abstractmethod
    def render( self, stream ):
        '''
        Render to the specified file object.
        '''
        pass

##############################################################################
class HtmlRunRenderer(RunRenderer):
    '''
    Render the results of run-time checking to an HTML document.
    '''
    def __init__( self, context, outParser, errParser ):
        self._context   = context
        self._outParser = outParser
        self._errParser = errParser

    def render( self, stream ):
        variables = {'context'   : self._context, \
                     'compiler'  : self._outParser.compiler,  \
                     'timestamp' : self._outParser.started, \
                     'events'    : self._errParser.getEvents()}

        environment = Environment(loader=PackageLoader('parserender', \
                                                       'templates'))
        template = environment.get_template( 'runlog.html' )
        print( template.render( variables ), file = stream )
