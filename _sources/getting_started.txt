.. getting_started

=========================
Installation Instructions
=========================

System requirements
-------------------

Recem [https://github.com/shoibalc/Recem/] has been tested against Python 2.7
on Linux (Ubuntu 11.10) and the latest
version of Coopr. Recem should work on all Linux computers but it might take
some work to get it to run on Windows or Mac OS X.

Installation Prerequisites
--------------------------

    * **Python 2.7**:

      Linux:
      Python should be available on every Linux/Unix system.
      Check the version using::

            $ python -V

      Install Python 2.7 from your Linux distribution repository if your
      version is older than that. On Ubuntu or Debian, the following command
      will do::

            $ sudo apt-get install python2.7

      Other Linuxes have similar commands: ``yum`` instead of ``apt-get`` on
      Fedora or RedHat.

      Windows:
      Install Python from this website:

            http://www.python.org/getit/windows/

      Macs:
      Macs come pre-installed with Python but the versions might be old.
      Install version 2.7 from this website:

            http://www.python.org/getit/mac/

    * **Coopr**:

      Coopr is A COmmon Optimization Python Repository - a collection of
      open source optimization packages in Python. Install Coopr from this
      website:

            https://software.sandia.gov/trac/coopr/wiki/GettingStarted

      We recommend that you use Options 2-4 (3 for Windows). Options
      2-4 set up a virtual Python environment which is a local and isolated
      Python environment that mimics a complete Python installation. We will
      install Recem in this virtual Python environment. The virtual Python
      environment comes with setuptools and pip - two very handy programs to
      install and manage other Python packages.

    * **COIN-OR Solvers**:

      COIN-OR is COmputational INfrastructure for Operations
      Research - an initiative to promote the development of open source
      software in operations research. We will use the solvers of COIN-OR.
      Install the solver executables from this website:

            http://www.coin-or.org/projects/CoinBinary.xml

      You should obtain the appropriate executables from the 'Download Binaries'
      link. Follow instruction their Installation help files to make sure that
      solvers are on your PATH (your operation system should know where to
      find them!).

Install Recem
--------------
We will install using the Python of the virtual Python installation of Coopr.
Assume that the Python executable of the virtual environment is
``/path/to/coopr/bin/python``.

Your's might be different. We assume that you have unzipped and saved the
Recem  zip file [https://github.com/shoibalc/Recem/zipball/master] in the
directory Recem. In this directory, type the following
commands::

     $ /path/to/coopr/bin/python setup.py develop

This will install a 'development' version of Recem. You can make changes,
uninstall Recem using::

     $ /path/to/coopr/bin/python setup.py develop --uninstall

and reinstall using::

     $ /path/to/coopr/bin/python setup.py develop

To see other options, type::

    $ /path/to/coopr/bin/python setup.py -h

Try it
-------

Type the following commands::

    $ cd recem/dice
    $ /path/to/coopr/bin/python runDICE1993.py
    $ /path/to/coopr/bin/python runDICE2007.py
