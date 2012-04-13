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


* MARKAL-TIMES type bottom-up model:

A MARKAL-TIMES bottom-up model is currently in development.
We expect to have a preliminary version for the USA and India.