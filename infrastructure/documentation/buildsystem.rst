The LFRic Build System
=======================

This is a quick introduction to the LFRic project's build system. Such a system
is intended to make developers' lives easier but to do so it requires a few
conventions to be followed.

.. note::
   The first LFRic application was and is developed in close collaboration with
   GungHo dynamical core developments, and was previously called `Dynamo`. As a
   result there are occasional references to `Dynamo` within the code base, the
   documentation and the Wiki pages. Such references will generally apply both
   to the development of the Dynamo application and to other applications
   (sometimes referred to as mini-apps) that use the LFRic infrastructure and
   LFRic build system.

.. note::
   The canonical version of this document exists as a reStructured text file
   in the repository at
   `source:LFRic/trunk/infrastructure/documentation/buildsystem.rst`:trac:.
   Branches which render it inaccurate should update the version of this
   document on that branch. The version displayed in the wiki is generated
   from the head of trunk.

Platform Identifiers
--------------------

In a number of places it is necessary to uniquely identify a target system.
This is achieved using a platform identifier.

Platform identifiers take the form <site>-<machine>. For example the Met
Office's XC40 is referred to as meto-xc40 while the current MONSooN machine
is known as monsoon-xc40. The previous MONSooN would have been monsoon-power7.
The workstation of someone working at Widget Corporation might be
widget-tigger, assuming their workstation was called Tigger.

Targets
-------

Everything is built via a top level ``make`` file. This offers a number of
options, referred to as targets, as described in the following sections.

Building
~~~~~~~~

To perform a full build simplt change to the top level of the working copy and
issue a ``make``::

  cd r1234_MyWorkingCopy
  make

Three build profiles are offered: ``full-debug``, ``fast-debug`` and
``production``. They are specified using the ``PROFILE`` variable passed to
make. For instance, to produce a build with all debugging enabled simply use
``make PROFILE=full-debug``. The default is ``fast-debug`` if none is
specified.

All three profiles insert debug symbols to the executable which provides
source file and line numbers in backtraces. They have no run-time impact other
than to increase the size of the executable.

``full-debug`` turns optimisation off and run-time checking, of e.g. array
out-of-bounds, on.

For those who need faster execution but also want bit reproducibility
``fast-debug`` limits optimisation to a safe level and turns off run-time
checking.

``production`` increases the optimisation level but does not guarantee
bit reproducibility. This is due to potential run-time optimisations,
particularly in the order of shared and distributed memory processing. We would
expect inter-run reproducibility but it can not be guaranteed.

Relocate Build Artifacts
^^^^^^^^^^^^^^^^^^^^^^^^

By default all build related files end up in the directory ``working`` within
the working copy. This is convenient for interactive development as it allows
the generated source to be examined. Unfortunately there are a number of
situations where this is less than ideal.

The database engine used to store dependency information interacts badly with
NFS exported mounts. If your working copy is stored on such a mount you may
find that builds are very slow. They will tend to stall during the construction
of ``programs.mk``.

Another filesystem which can have problems is Lustre, often used on
supercomputers. It is tuned for large parallel access to a few files, not
the many serial accesses to small files involved in compiling. Thus a working
copy held on such a file system will be slow to compile and place unreasonable
demands on the file server.

Both these problems can be aleviated by relocating the build directory.

This is achieved by simply exporting the ``WORKING_DIR`` environment
variable. It should contain an absolute path to a suitible temporary space in
which to perform the compile

.. note::
  Regardless of where the working directory is the finished build products will
  end up in directories ``bin`` and ``tests`` in the root of the working copy.

Verbose Building
^^^^^^^^^^^^^^^^

By default the build system will suppress as much output as possible to reduce
clutter. When resolving problems it may be useful to see the output. If you want
to see the actual commands being run by the build set the ``VERBOSE`` variable::

  make VERBOSE=1

MPI
^^^

If you are attempting to build with MPI support do *not* set ``FC`` to
``mpif90``. For a start some compilers (Cray, IBM) deal with MPI themselves.

Instead define the variable ``LDMPI`` with the command to link against the MPI
libraries. This will generally be ``mpif90``.

Linking
^^^^^^^

By default the build system will link dynamically. If you want a statically
linked binary then you should pass the ``LINK_TYPE`` variable::

  make LINK_TYPE=static

Static linking is the default on Cray systems as they do not seem to play well
with dynamic linking.

Cleaning
^^^^^^^^

As with many build systems it is possible to ``make clean`` to delete all build
artefacts. These include working files and complete executable binaries.

Testing
^^^^^^^

Unit tests will be built and run as part of a normal build, so there is no
need to worry about that. The test suite, on the other hand, must be
manually invoked.

When ``make test-suite`` is used a number of instances of Rose Stem will be
launched. The environment variable ``TEST_SUITE_TARGETS`` holds a space
separated list of the platform identifiers. These identify the targets to be
used for test suite runs.

Configurations are held in ``rose-stem/opt``. Each filename has the form
``rose-suite-<platform id>.conf`` using the platform identifiers described
above.

For those using a Met Office module collection the core module will set this up
for you. e.g. On the desktop do::

  module load common-environment/lfric

For further information an testing see `Dynamo/Testing`:trac:.


Choosing a Compiler
-------------------

A number of compilers are supported by the build system. To choose between them
simply set the ``FC`` environment variable. e.g. ``export FC=ifort`` This is
often done by an automated system, as it is for users of the Met Office LFRic
module system.

Each compiler needs a different set of arguments. These are defined by the
LFRic build system in ``infrastructure/build/fortran/<compiler>.mk``. Each
file is named after the command used to launch the compiler.

A number of universal compiler options and build system variables are defined
first.

=========================  =================================================================================
Variable                   Purpose
=========================  =================================================================================
``<compiler>_VERSION``     Holds the version of the running compiler. This is held as a zero padded integer. e.g. 4.9.2 becomes 040902
``F_MOD_DESTINATION_ARG``  The argument to be used to tell the compiler where to put ``.mod`` files.
``OPENMP_ARG``             The argument passed to enable interpretation of OpenMP directives.
``DEPRULE_FLAGS``          Additional arguments to pass to the dependency generator script.
=========================  =================================================================================

Notice that the ``FFLAGS`` variable is not modified directly. Instead a number
of argument groups are defined.

=============================  ===========================================================================================================
Argument Group                 Purpose
=============================  ===========================================================================================================
``FFLAGS_COMPILER``            Control the way the compiler works. For instance mdoifying unexpected behavior.
``FFLAGS_NO_OPTIMISATION``     Turn off all optimisation.
``FFLAGS_SAFE_OPTIMISATION``   Enable only bit reproducibility safe optimisation.
``FFLAGS_RISKY_OPTIMISATION``  Enable optimisations which may affect bit reproducibility.
``FFLAGS_DEBUG``               Enable the generation of debugging symbols, tracebacks and such.
``FFLAGS_WARNINGS``            Control warning handling. For development we like to receive all warnings and have them act as errors.
``FFLAGS_INIT``                Turn on any declaration initialisation the compiler may support to fill declared memory with a known value.
``FFLAGS_RUNTIME``             Turn on any run-time checking the compiler may support. Array bound checking for instance.
=============================  ===========================================================================================================

The same is done with ``LDFLAGS``.

=========================  =======================================================================================================
Variable                   Purpose
=========================  =======================================================================================================
``LDFLAGS_COMPILER``       Some linkers need to be told to retain debug symbols. That could be done here.
``LD_COMPILER_LIBRARIES``  Specify any special libraries needed at link stage as a space separated list.
                           These are not filenames, drop the "lib" from the front.
=========================  =======================================================================================================


Automatic Dependency Analysis
-----------------------------

A knowledge of dependencies between source files allows only those things
affected by a change to be recompiled. Maintaining such a dependency list by
hand is tedious and error prone. This is why the build system will generate
and maintain dependency information automatically.

The dependency analyser scans Fortran source looking for "use" statements
to determine the prerequisites of a source file. It relies on source files
having the same name as the module they contain, which implies that they contain
only one module.

e.g. if the module is called "my_special_stuff_mod" then the file name
should be "my_special_stuff_mod.f90". If this is not the case these modules will
be rebuilt every time.

The current dependency analyser works well for most common uses, but it will
fail when you rename files. If you rename a file, run ``make clean`` prior
to rebuilding. This will force the dependency database to be rebuilt.


External Libraries
------------------

If you want to make use of an external library then some simple edits to
``Makefile`` will be required.

There are two critical variables: ``IGNORE_DEPENDENCIES`` and
``EXTERNAL_*_LIBRARIES``.

``IGNORE_DEPENDENCIES`` is a space-separated list of modules which should be
ignored by the dependency analyser when discovered in ``use`` statements. This
prevents the dependency analysis trying to rebuild your library. For example,
the infrastructure makes use of the ESMF library: the name to add would be
"esmf" since we ``use esmf``.

``EXTERNAL_*_LIBRARIES`` is a space-separated list of library names to be passed
to the linker using "little l" arguments. Libraries which are available as a
".so" file should be listed in ``EXTERNAL_DYNAMIC_LIBRARIES`` while those only
available as a ".a" file are listed in ``EXTERNAL_STATIC_LIBRARIES``. Returning
to the ESMF example, the library file is ``libesmf.so`` so the string ``esmf``
would be added to ``EXTERNAL_DYNAMIC_LIBRARIES``.

In addition to adding knowledge of the library to the build system you have to
make it findable. If you are using a module from the Met Office module
collection then this is handled for you automatically. If you are not then you
will need to do a little additional work.

The path to the directory containing the ``mod`` files should be added to
``FFLAGS`` in "big I" notation. e.g. ``export FFLAGS="$FFLAGS -I/path/to/mods"``

The path to the directory containing the library files should be added to
``LDFLAGS`` in "big L" notation.
e.g. ``export LDFLAGS="$LDFLAGS -L/path/to/libs"``

This will get you compiling. To run the resulting executable you have to make
sure the run-time linker can find the library files. To achieve this modify
``LD_LIBRARY_PATH``.
e.g. ``export LD_LIBRARY_PATH=/path/to/libs:$LD_LIBRARY_PATH``


Automatically Generated Code
----------------------------

The LFRic infrastrucuture make use of the PSyclone tool, created as part of the
Gung Ho project, to automatically generate some of the code. PSclone is
designed to generate the Parallel Systems layer (PSy-layer) code, but it also
modifies the algorithm code as written by science developers.

Algorithm routines in ``gungho/source/algorithm`` are written with a ``.x90``
extension. From these both compilable Fortran source for the algorithm and
corresponding PSy-layer are generated. Normally this process is handled for
you by the build system.

Where PSyclone does not yet support a feature that you wish to use, you can
override the automatically generated code with a manually written PSy-layer
module.

The simplest way to achieve this is to start by building Gung Ho with a stock
checkout. This will generate all the PSy-layer modules as normal. Then move the
one you wish to override from "working/gungho/algorithm" to
``gungho/source/psy`` and modify it according to your needs. The file will have
the form ``*_psy.f90``. Algorithm files are automatically generated in all cases
so should not be moved in this way.

.. attention::
   You should ensure that the LFRic and PSyclone developers are aware of your
   need to modify the PSy-layer code, to ensure that your changes fit in with
   the ongoing development of LFRic and PSyclone.

Optimisation
~~~~~~~~~~~~

PSyclone is able to make use of scripts to apply optimisations to the code as
it is generating it.

It is necessary for PSyclone to know which platform you are building on in
order to select the correct optimisation scripts. This is achieved through the
`OPTIMISATION_PROFILE` environment variable. It should contain a single
platform identifier following the convention outlined above.

The build system will look for ``optimisation/<platform id>/global.py`` which
will be applied to all algorithm files. Platform identifiers are specified
above.

If a file ``optimisation/<platform id>/<algorithm>.py`` exists it will be used
in preference to the global script. The algorithm name is taken from
``gungho/source/algorithms/<algorithm>.x90``.

UM physics codes
~~~~~~~~~~~~~~~~

We attempt to develop LFric with single source physics taken from the UM
repository.  In order to build in the UM code, fcm make is used to extract
and preprocess the code; this is then rsync'd with the working build directory
such that the LFRic build system can proceed with analysing and building this 
code.  Since this process results in the build system analysing a considerable
amount of additional code, the variable  ``UM_PHYSICS`` is passed to the make
invocation to indicate this process is wanted, e.g.

  ``make UM_PHYSICS=1 build-gungho``

The make procedure will then carry out the ``fcm make`` invocation and then 
subsequently rsync the extracted and preprocessed code to the ``working`` directory 
tree. It is the intention that the UM code for the build is kept separate from the main 
LFRic source and any modifications on the UM side should be made through the branches 
incorporated at the fcm make stage (these could be a separate working copy). To change the 
UM branches incorporated into the build, modify the ``um_sources`` environment variable in the 
``set_environment.sh`` file.

Since the UM code will continue to evolve and we will want to source difference 
versions/branches from the UM repository, the environment variables needed by the fcm make
command (including the paths to the repository/branches/working) are set in 
``um_physics/set_environment.sh``.  Typically, this script will be maintained in the 
LFRic repository in but can be overridden if a different UM source/configuration is required.

fcm-make
~~~~~~~~

This section has been rendered inacurate and the associated apability
inoperative by the replacement of the build system. It is left here for
reference until such times as it can be redeveloped.

The Dynamo code, together with PSyclone and pFUnit, can be extracted using 
fcm-make. This can extract direct to the local machine without the need for a
mirroring step.

  ``fcm make -f fcm-make/<site-platform-compiler>/<optimisation>.cfg``

This will extract the code to a directory named "extract". Because fcm-make
extracts code in a different directory structure than the compilation system
expects, the following environment variables need to be set before invoking
``make``. Note that the final three environment variables are only required to 
run the unit tests.

=============================  ===========================================================================================================
Environment Variable           Setting
=============================  ===========================================================================================================
``DYNAMO_BUILD_ROOT``          $PWD/extract/dynamo
``PSYCLONE``                   python2.7 $PWD/extract/psyclone/src/generator.py
``PSYCLONE_DIR``               $PWD/extract/psyclone
``PYTHONPATH``                 $PWD/extract/psyclone/f2py_93:$PWD/extract/psyclone/src:$PYTHONPATH
``PFUNIT_SOURCE_DIR``          $PWD/extract/pfunit
``PFUNIT_BUILD_DIR``           $PWD/extract/pfunit
``INSTALL_DIR``                $PWD/extract/dynamo/pfunit-install
``F90``                        ifort
``F90_VENDOR``                 Intel
=============================  ===========================================================================================================

