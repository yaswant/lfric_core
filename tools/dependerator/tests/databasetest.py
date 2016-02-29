#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

import unittest
import os.path
import shutil
import sqlite3
import tempfile

import dependencyanalysis.database

###############################################################################
class DatabaseTest( unittest.TestCase):
    ###########################################################################
    def setUp( self ):
        self._scratchDirectory = tempfile.mkdtemp()
        self._dbFilename = os.path.join( self._scratchDirectory, 'test.db' )
        self._database = sqlite3.connect( self._dbFilename )

    ###########################################################################
    def tearDown( self ):
        self._database.close()
        shutil.rmtree( self._scratchDirectory )

    ###########################################################################
    def testConstructor( self ):
        # First we test an empty DB
        uut = dependencyanalysis.database.Dependencies( self._dbFilename )

        cursor = self._database.cursor()

        cursor.execute( 'SELECT * FROM provides' )
        row = cursor.fetchall()
        self.assertEqual( len(row), 0 )

        cursor.execute( 'SELECT * FROM dependencies' )
        row = cursor.fetchall()
        self.assertEqual( len(row), 0 )

        cursor.execute( 'SELECT * FROM programs' )
        row = cursor.fetchall()
        self.assertEqual( len(row), 0 )

        # Now that we've create a DB make sure we don't nuke it.
        # First put some dummy data in.
        cursor.execute( 'INSERT INTO provides VALUES ( "foo.x", "fu" )' )
        self.assertEqual( cursor.rowcount, 1 )
        cursor.execute( 'INSERT INTO dependencies VALUES ( "bar", "qux" )' )
        self.assertEqual( cursor.rowcount, 1 )
        cursor.execute( 'INSERT INTO programs VALUES ( "baz" )' )
        self.assertEqual( cursor.rowcount, 1 )

        # Then create a new object and check the DB still contains what we
        # think it should
        uut = dependencyanalysis.database.Dependencies( self._dbFilename )

        cursor.execute( 'SELECT * FROM provides' )
        row = cursor.fetchall()
        self.assertEqual( row, [('foo.x', 'fu' )] )

        cursor.execute( 'SELECT * FROM dependencies' )
        row = cursor.fetchall()
        self.assertEqual( row, [('bar', 'qux')] )

        cursor.execute( 'SELECT * FROM programs' )
        row = cursor.fetchall()
        self.assertEqual( row, [('baz',)] )

    ###########################################################################
    def testAddProgram( self ):
        uut = dependencyanalysis.database.Dependencies( self._dbFilename )
        uut.addProgram( 'foo', 'bar.x' )

        row = self._database.execute( 'select * from provides' ).fetchall()
        self.assertEqual( row, [('bar.x', 'foo')] )

        row = self._database.execute( 'select * from programs' ).fetchall()
        self.assertEqual( row, [('foo',)] )

    ###########################################################################
    def testAddModule( self ):
        uut = dependencyanalysis.database.Dependencies( self._dbFilename )
        uut.addModule( 'foo', 'bar.x' )

        row = self._database.execute( 'select * from provides' ).fetchall()
        self.assertEqual( row, [('bar.x', 'foo')] )

    ###########################################################################
    def testAddDependency( self ):
        uut = dependencyanalysis.database.Dependencies( self._dbFilename )
        uut.addDependency( 'foo', 'bar' )

        row = self._database.execute( 'select * from dependencies' ).fetchall()
        self.assertEqual( row, [('foo', 'bar')] )

        # Add a second dependency on the same dependor
        uut.addDependency( 'foo', 'qux' )

        row = self._database.execute( 'select * from dependencies' ).fetchall()
        self.assertEqual( row, [('foo', 'bar'), ('foo', 'qux')] )

    ###########################################################################
    def testRemoveSourceFile( self ):
        uut = dependencyanalysis.database.Dependencies( self._dbFilename )
        self._populateDB()

        # Delete something from the middle of a dependency chain
        uut.removeSourceFile( 'bar.x' )

        rows = self._database.execute( 'SELECT * FROM provides' ).fetchall()
        self.assertEqual( rows, [('foo.x', 'foo'), \
                                 ('baz.x', 'baz'), \
                                 ("qux.x", "qux"), \
                                 ("fred.x", "fred"), \
                                 ("wilma.x", "wilma")] )

        rows = self._database.execute( 'SELECT * FROM dependencies' ).fetchall()
        self.assertEqual( rows,[("foo", "bar"), \
                                ("foo", "baz"), \
                                ("fred", "wilma")] )

        rows = self._database.execute( 'SELECT * FROM programs' ).fetchall()
        self.assertEqual( rows, [('foo',), ("fred",)] )

    ###########################################################################
    def testGetDependencySources( self ):
        uut = dependencyanalysis.database.Dependencies( self._dbFilename )
        self._populateDB()

        dependencies = [tuple(row) for row in uut.getDependencySources( 'foo.x' )]
        self.assertEqual( len(dependencies), 2 )
        self.assertItemsEqual( dependencies, [('foo.x', 'bar.x', 'bar'), \
                                              ('foo.x', 'baz.x', 'baz')] )

    ###########################################################################
    def testProgramSourceList( self ):
        uut = dependencyanalysis.database.Dependencies( self._dbFilename )
        self._populateDB()

        sources = [tuple(row) for row in uut.getProgramSourceList()]
        self.assertEqual( len(sources), 2 )
        programs = []
        expected = { 'foo':['foo.x', 'bar.x', 'baz.x', 'qux.x'], \
                     'fred':['fred.x', 'wilma.x'] }
        for program, files in sources:
            programs.append( program )
            self.assertIn( program, expected.keys() )
            self.assertItemsEqual( expected[program], files )
        self.assertItemsEqual( ['foo', 'fred'], programs )

    ###########################################################################
    def testProgramSources( self ):
        uut = dependencyanalysis.database.Dependencies( self._dbFilename )
        self._populateDB()

        sources = uut.getProgramSources()
        self.assertEqual( {'foo':'foo.x', 'fred':'fred.x'}, sources )

    ###########################################################################
    def testGetAllFileDependencies( self ):
        uut = dependencyanalysis.database.Dependencies( self._dbFilename )
        self._populateDB()

        found = []
        expected = { 'foo.x':['bar.x', 'baz.x'], \
                     'bar.x':['qux.x'], \
                     'fred.x':['wilma.x'] }
        for dependor, dependees in uut.getAllFileDependencies():
            found.append(dependor)
            self.assertIn( dependor, expected.keys() )
            self.assertItemsEqual( expected[dependor], dependees )
        self.assertItemsEqual( ['foo.x', 'bar.x', 'fred.x'], found )

    ###########################################################################
    def _populateDB( self ):
        affected = self._database.executemany( 'INSERT INTO provides VALUES (?,?)', \
                                               [("foo.x", "foo"), \
                                                ("bar.x", "bar"), \
                                                ("baz.x", "baz"), \
                                                ("qux.x", "qux"), \
                                                ("fred.x", "fred"), \
                                                ("wilma.x", "wilma")] ).rowcount
        self.assertEqual( affected, 6 )

        affected = self._database.executemany( 'INSERT INTO dependencies VALUES (?,?)', \
                                               [("foo", "bar"), \
                                                ("foo", "baz"), \
                                                ("bar", "qux"), \
                                                ("fred", "wilma")] ).rowcount
        self.assertEqual( affected, 4 )

        affected = self._database.executemany( 'INSERT INTO programs VALUES (?)', \
                                               [("foo",), ("fred",)] ).rowcount
        self.assertEqual( affected, 2 )

        self._database.commit()

if __name__ == '__main__':
    unittest.main()
