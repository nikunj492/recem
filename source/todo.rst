.. todo

================
Project Contents
================

* DICE:

Two versions of `Wiliam Nordhaus <http://nordhaus.econ.yale.edu/>`_' DICE
integrated assessment model - DICE1993 and DICE2007 have been coded in Pyomo
as proof of concept examples. The dice models can be created by importing
the recem.dice package in a python program::

    from coopr.pyomo import *
    from recem.dice import createDICE1993, createDICE2007

    dice1 = createDICE1993()
    dice2 = createDICE2007()

The python scripts **runDICE1993.py** and **runDICE2007.py** provide example
 programs for their use.
::

    $ python runDICE1993.py

will generate the file dice1993.txt, which contains the results of the runs of
BASE, MARKET, OPT_CONT and CONCENT scenarios of DICE1993. The file is a csv file with
'space' or [_] as the delimiter, any spreadsheet program should be able to open it.
Similarly::

    $ python runDICE2007.py

will generate the file dice2007_optimum.txt, which contains the results for the
OPTIMUM scenario of DICE2007.

* MARKAL-TIMES type bottom-up model:

A MARKAL-TIMES bottom-up model is currently in development.
We expect to have a preliminary version for the USA and India.