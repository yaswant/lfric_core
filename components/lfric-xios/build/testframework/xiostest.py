#!/usr/bin/env python3
##############################################################################
# (C) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

from __future__ import print_function

from __future__ import absolute_import
from abc import ABCMeta, abstractmethod
import os
from pathlib import Path
import sys
from typing import List
import six
from testframework import MpiTest


##############################################################################
class LFRicXiosTest(six.with_metaclass(ABCMeta, MpiTest)):
    """
    Base for LFRic-XIOS integration tests.
    """

    def __init__(self, command=sys.argv[1], processes=1):
        """ """
        super(LFRicXiosTest, self).__init__(command, processes)
        self.xios_out: List[XiosOutput] = []
        self.xios_err: List[XiosOutput] = []

    def post_execution(self, return_code):
        """
        Cache XIOS logging output for analysis.
        """

        for proc in range(self._processes):
            self.xios_out.append(XiosOutput(f"xios_client_{proc}.out"))
            self.xios_err.append(XiosOutput(f"xios_client_{proc}.err"))


class XiosOutput:
    """
    Simple class to hold XIOS output log information
    """

    def __init__(self, filename):
        self.path: Path = Path(os.getcwd()) / Path(filename)

        with open(self.path, "rt") as handle:
            self.contents = handle.read()

    def exists(self):
        """
        Checks if log output file exists
        """
        return os.path.exists(self.path)
