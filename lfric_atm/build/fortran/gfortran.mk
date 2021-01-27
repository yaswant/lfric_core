##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Various things specific to the GNU Fortran compiler.
##############################################################################

$(info Project specials for GNU compiler)

export FFLAGS_UM_PHYSICS = -fdefault-real-8
# lfric_atm dependencies contain code with implicit lossy conversions and
# unused variables.
# We reset the FFLAGS_WARNINGS variable here in order to prevent
# -Werror induced build failures.
FFLAGS_WARNINGS          = -Wall -Werror=character-truncation -Werror=unused-value
