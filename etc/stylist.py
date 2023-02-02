##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# Configures the code style used for automatic checking using Stylist.
#
# Stylist is our code style checking utility. It may be found at
# <https://github.com/MetOffice/stylist> with documentation at
# <https://metoffice.github.io/stylist>.
#
# For the full project code style check this page:
# <https://code.metoffice.gov.uk/trac/lfric/wiki/LFRicTechnical/CodingStandards>
#
from re import compile as re_compile

from stylist.fortran import (
    FortranCharacterset,
    MissingImplicit,
    IntrinsicModule
)
from stylist.rule import TrailingWhitespace
from stylist.source import (
    FilePipe,
    FortranPreProcessor,
    FortranSource,
    PFUnitProcessor
)
from stylist.style import Style


if __name__ == '__main__':
    raise Exception("This is a Stylist configuration file, not to be run")

# Define the rules which make up our style...
#
infrastructure = Style(
    TrailingWhitespace(),
    FortranCharacterset(),
    MissingImplicit(),
    IntrinsicModule()
)

# Define additional file type processing pipelines
#
pf = FilePipe(FortranSource, PFUnitProcessor, FortranPreProcessor)
