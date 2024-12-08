.. |ech2o| replace:: EcH\ :sub:`2`\ O

Installation
=============

The following page contains information on how to compile the latest version of |ech2o|-iso in Linux, Windows and Mac, whose executable is ``ech2o_iso``.
Instructions to install pre-compiled Linux and Windows version will be available soon.
For now, pre-compiled executables are provided in the ``CaseStudy`` folder ; note however that the Linux and Mac executable are platform-dependent and might not run out of the box!

|ech2o|-iso uses PCRASTER as data pre- and post-processor. Please install PCRASTER free of charge from `here <https://pcraster.geo.uu.nl/downloads/latest-release/>`_.

1. Compilation instructions for |ech2o|-iso using gcc
-----------------------------------------------------
    
 The current version of |ech2o|-iso does not have a configure script. The Makefile has been generated for the gnu c++ compiler and does not check for dependencies. MINGW is necessary to compile in Windows. 


1. Dependencies
^^^^^^^^^^^^^^^^

* Change directory to your workspace and clone the latest version of the source files from the git repository:

::

   $ git clone https://bitbucket.org/scicirc/ech2o-iso.git

* Install the armadillo development files (version 7 or higher), either compiling and installing from source (`here <http://arma.sourceforge.net/download.html>`_) or from the package manager of your Linux distribution (or from homebrew on Mac).

* Install the boost development files, either compiling and installing from source (`here <http://www.boost.org>`_) or from the package manager of your Linux distribution. For Mac users, note that the boost libraries installed with homebrew are compiled with Apple Clang (which should not be used here, as it does not handle OpenMP), and will not link at the end |ech2o|-iso compiling if the latter is done with gcc. One workaround is either to compile boost from source using gcc (see discussion `here <https://stackoverflow.com/questions/25346443/how-to-install-boost-with-specified-compiler-say-gcc>`_) or compile boost and |ech2o|-iso using LLVM Clang.

* Precompiled versions of the libcsf dependency for Linux, Windows and Mac are included in the ``lib`` folder. The compilation was carried assuming little endian 64 bit architectures.

  If the linker complains, the library may need to be compiled for your system. Please, clone the source code from 
    
::
   
   $ git clone https://bitbucket.org/maneta/libcsf.git
   

or download from
   
::
   
   $  https://sourceforge.net/p/pcraster/rasterformat/ci/master/tree/
   
and compile from source. Then replace the old libcsf64 library in the ``lib`` directory with the newly compiled library. Make sure you change the name of the new library so it has the same name as the old one. 
   

2. Making ``ech2o_iso``
^^^^^^^^^^^^^^^^^^^^^^^^

* Open a command line terminal

* Change to the ``Release-*`` folder within the source folder, where ``*`` is your OS type: Linux, Mac or Windows.

* Type ``make all`` (or ``make ech2o_iso``) to make the source. Use ``make clean`` to reset (it removes the *.o, *.d and the ``ech2o_iso`` executable).

.. Note::
   If compiling for Windows, edit the objects.mk file and substitute item ``-lcsf64`` for ``-llibcsf64`` so that ``make`` will link against the correct static library. Save and close the editor
   If compiling for Mac, edit the objects.mk file and substitute item ``-lcsf64`` for ``-lcsfosx`` so that ``make`` will link against the correct static library. Save and close the editor

3. Making ``asc2c``
^^^^^^^^^^^^^^^^^^^^

* Open a command line terminal 
 
* Change directory to your workspace and clone the latest version of the source files from the git repository:

::

   $ git clone https://bitbucket.org/maneta/asc2c.git

* Change directory into the source folder and type ``make`` to make the code. 

  
  
Contact
----------

If you need assistance compiling the source, contact sylvain[dot]kuppel[at]ird[dot]fr or marco[dot]maneta[at]umontana[dot]edu 

If you find this documentation to be incomplete, please file a ticket in the appropriate issue tracker:

* ech2o_iso compilation issues:  https://bitbucket.org/scicirc/ech2o-iso/issues
* asc2c compilation issues:  https://bitbucket.org/maneta/asc2c/issues
  
