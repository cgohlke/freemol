FreeMOL Design Document (draft ideas subject to further refinement)

Goal: 

  Aggregate a managable set of open-source molecular software tools
  together into a maintainable form which can support integration with
  PyMOL and other interactive applications via quality-controlled
  platform-specific binary distributions.

Structure:

  "src" top-level directory with distinct subdirectories containing
  configure scripts and source code for the various FreeMOL packages

  "freemol" top-level directory is a prototype framework for a freemol
  binary distribution that would contain executables and other
  required runtime files in a typical "unix-/usr-like" tree:
  (freemol/bin, freemol/lib, freemol/share, freemol/man, and so
  forth).  The FREEMOL environment variable should be set to point at
  this directory.  Please note that platform-specific binary
  executables should not be checked into subversion into this
  prototype directory tree.

Compilation:

  Assuming the the FREEMOL environment variable has been defined,
  each package can be compiled and installed into $FREEMOL by cd'ing
  to its souce directory and issuing the following sequence of commands:

    ./configure
    make
    make install
    
  Followed by a 

    make clean

  to get rid of object code, intermediate files, etc.

Deployment:

  A FreeMOL binary distribution would be a tar or zipped archive of
  the $FREEMOL folder containing executables compiled for a specific
  architecture.  

  Installation will amount to extraction of the archive and setting
  of the FREEMOL environment variable to that PyMOL and other packages
  can subquently locate and rely upon the freemol executables. 

Pragmatic decisions:

  Rely upon unix-like shell-based compilation (MinGW/Cygwin on Win32)

  Stick with standalone, command-line only packages which have few if
  any dependencies

Confirmed FreeMOL Packages:

  mengine 
  mpeg_encode 
  apbs 
  pdb2pqr
  
Proposed FreeMOL Packages

  python
  rdkit
  openbabel








