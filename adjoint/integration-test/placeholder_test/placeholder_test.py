#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Run a placeholder integration test
'''

import os
import re
import sys

from testframework import LFRicLoggingTest, TestEngine, TestFailed


class PlaceHolderTest(LFRicLoggingTest):
    '''
    Run a placeholder integration test
    '''
    _RESULT_PATTERN = re.compile( r'\bResult\s*=\s*([-0-9.]+\b)', re.IGNORECASE)
    def __init__(self):
        if 'MPIEXEC_BROKEN' in os.environ:
            PlaceHolderTest.set_mpiexec_broken()
        super(PlaceHolderTest, self).__init__([sys.argv[1],
                                              'placeholder_test_configuration.nml'],
                                              processes=1,
                                              name='PlaceholderTest.Log')

    def test(self, return_code, out, err):
        '''
        Error messages if the test failed to run
        '''
        if return_code != 0:
            message = 'Test program failed with exit code: {code}'
            raise TestFailed(message.format(code=return_code),
                             stdout=out, stderr=err,
                             log=self.getLFRicLoggingLog())

        for line in out.split("\n"):
            match = PlaceHolderTest._RESULT_PATTERN.search(line)
            if match:
                expected = 0.0
                if (float(match.group(1)) - expected) < 0.001:
                    return 'Science and tech worked together'
                else:
                    raise TestFailed('Expected >{}< but found >{}<' \
                                     .format(expected, match.group(1)))

        # Case for when there is no match
        raise TestFailed('Unable to find result in output: ' + out)

if __name__ == '__main__':
    TestEngine.run( PlaceHolderTest() )
