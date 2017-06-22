#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
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
    def render( self, fileList ):
        '''
        Render to the specified file object.
        '''
        return ''

##############################################################################
class HtmlIndexRenderer(IndexRenderer):
  def __init__( self ):
    pass

  def render( self, title, fileList ):
    variables = {'title' : title, \
                 'files' : fileList}
    environment = Environment( loader=PackageLoader( 'parserender', \
                                                     'templates' ) )
    template = environment.get_template( 'simpleindex.html' )
    return template.render( variables )

##############################################################################
class SiteIndexRenderer(object):
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
class HtmlSiteIndexRenderer(SiteIndexRenderer):
    '''
    Render the index page to an HTML document.
    '''
    def __init__( self, indexer, suiteUrl ):
        self._indexer  = indexer
        self._suiteUrl = suiteUrl

        if suiteUrl is not None:
          if self._suiteUrl.endswith('/'):
            self._suiteUrl = self._suiteURL[:-1]

    def render( self, stream ):
        variables = {'cronout'       : self._indexer.cronOut,
                     'crontimestamp' : self._indexer.cronTimestamp,
                     'oldthreshold'  : datetime.datetime.utcnow()
                                       - datetime.timedelta( hours=48 ),
                     'suiteurl'      : self._suiteUrl,
                     'tree'          : self._indexer.getTree()}

        environment = Environment( extensions=['jinja2.ext.do',
                                               'jinja2.ext.loopcontrols'],
                                   loader=PackageLoader( 'parserender',
                                                         'templates' ) )
        template = environment.get_template( 'siteindex.html.jinja' )
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
    def __init__( self, context, statusParser, outParser, errParser ):
        self._context      = context
        self._statusParser = statusParser
        self._outParser    = outParser
        self._errParser    = errParser

    def render( self, ignoreCodes, stream ):
        if ignoreCodes:
            eventBuckets = self._errParser.getEvents( ignoreCodes=ignoreCodes, \
                                                      groupBy='filename' )
        else:
            eventBuckets = self._errParser.getEvents( groupBy='filename' )

        variables = {'context'      : self._context, \
                     'compiler'     : self._outParser.compiler,  \
                     'timestamp'    : self._statusParser.started, \
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
