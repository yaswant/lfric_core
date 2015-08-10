#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
'''!
Library of tools for handling Subversion repositories and working copies.
'''

from __future__ import print_function

import calendar
import datetime
import subprocess
import xml.sax
import xml.sax.handler

__version__ = '1.0.0'

def subversionInfo(repoURL):
    '''!
    Get repository information.

    The returned dictionary contains the following:

| Key                  | Value                                                 |
| -------------------- | ----------------------------------------------------- |
| latest-revision      | The revision which most recently modified the URL.    |
| latest-revision-time | The (Unix) timestamp of the most recent modification. |


    @param repuURL (str) The URL of the repository.
    @return (dict) Map of parameter to value.
    '''
    infoHandler = _SubversionInfoHandler()
    process = subprocess.Popen( ['svn', 'info', '--xml', repoURL], \
                                stdout=subprocess.PIPE )
    xml.sax.parse(process.stdout, infoHandler)

    info = { 'latest-revision'      : infoHandler.revision, \
             'latest-revision-time' : infoHandler.timestamp }
    return info

def subversionLog(location):
    '''!
    Get the change log for a particular directory.

    @param location (str) Either a repository URL or working copy path.
    @return (dict) Map of timestamp (Unix time) to log message.
    '''
    logHandler = _SubversionLogHandler()
    process = subprocess.Popen( ['svn', 'log', '--xml', location], \
                                stdout=subprocess.PIPE )
    xml.sax.parse( process.stdout, logHandler )
    return logHandler.events

class SubversionException(Exception):
    '''!
    Thrown by this module to indicate fault conditions.
    '''
    pass

def _parseTimestamp(stampString):
    '''!
    Convert a Subversion timestamp into Unix time.

    @param stampString (str) The timestamp as culled from the XML.
    @return (int) Number of seconds since Unix epoch.
    '''
    timestamp = datetime.datetime.strptime(stampString, \
                                           '%Y-%m-%dT%H:%M:%S.%fZ')
    return calendar.timegm(timestamp.utctimetuple())

class _XMLHandler( xml.sax.handler.ContentHandler, object ):
    '''!
    Provides basic handling of SAX events. It is intended that this class will
    be specialised with your own specifics but that you will call down to it
    for the basics.
    '''
    def __init__( self ):
        '''!
        Initialiser.
        '''
        self._tagStack = []

    def documentTag( self ):
        '''!
        Get the root tag of the document. Override this for a specific
        document.

        @return (str) The name of the document tag.
        '''
        raise NotImplementedError('Interface not fully implemented')

    def startElement( self, name, attrs ):
        '''!
        Track and sanity check the tag heirarchy.
        '''
        self._tagStack.append(name)
        self._content = ''

        if ( len ( self._tagStack ) == 1 ) and ( name != self.documentTag() ):
            raise SubversionException( \
                'Parsing {} document but found {}'.format( self.documentTag(), \
                                                           name ) )

    def endElement( self, name ):
        '''!
        Track and sanity check the tag heirarchy.
        '''
        expected = self._tagStack.pop()
        if name != expected:
            raise SubversionException( \
                'Found end of {} but expected {}'.format( name, expected ) )

    def characters( self, content ):
        '''!
        Concatenate tag content into a single string.
        '''
        self._content += content

class _SubversionInfoHandler( _XMLHandler ):
    '''!
    Handle SAX events coming from parsing the output of "svn info".

    Upon completion the object will have two public members: revision and
    timestemp. These hold the revision number of the most recent commit and
    time time it was made.
    '''
    def documentTag( self ):
        return 'info'

    def startElement( self, name, attrs ):
        super( _SubversionInfoHandler, self).startElement( name, attrs )

        if name == 'commit':
            self.revision = attrs[ 'revision' ]

    def endElement( self, name ):
        super( _SubversionInfoHandler, self).endElement( name )

        if name == 'date':
            if self._tagStack[ -1 ] == 'commit':
                timestamp = datetime.datetime.strptime(self._content, \
                                                    '%Y-%m-%dT%H:%M:%S.%fZ')
                self.timestamp = calendar.timegm(timestamp.utctimetuple())

class _SubversionLogHandler( _XMLHandler ):
    '''!
    Handle SAX events coming from parsing the output of "svn log".

    Upon completion the object will have a public member "events" which is a
    map of revision timestamp to commit message.
    '''
    def __init__(self):
        super( _SubversionLogHandler, self).__init__()

        self.events = {}

    def documentTag( self ):
        return 'log'

    def startElement( self, name, attrs ):
        super( _SubversionLogHandler, self).startElement( name, attrs )

        if name == 'logentry':
            self._revision = attrs[ 'revision' ]

    def endElement( self, name ):
        super( _SubversionLogHandler, self).endElement( name )

        if name == 'date':
            timestamp = datetime.datetime.strptime(self._content, \
                                                   '%Y-%m-%dT%H:%M:%S.%fZ')
            self._timestamp = calendar.timegm(timestamp.utctimetuple())
        elif name == 'msg':
            self._message = self._content.strip()
        elif name == 'logentry':
            self.events[ self._timestamp ] = self._message

if __name__ == '__main__':
    pass
