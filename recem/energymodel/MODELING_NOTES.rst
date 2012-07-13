#==============================================================================
# RECEM: Repository of Energy, Climate, and Economics Models
#
# Copyright (c) 2012 Princeton University
# Author: Shoibal Chakravarty
# This software is distributed under the BSD License.
# You may obtain a copy of the License at
# http://www.opensource.org/licenses/bsd-license.php
#==============================================================================


OPTIMIZATION MODELS:
====================

In optimization models a decision making problem is solved. The decision is
usually about what investments to make or resources to run in order to meet a
constraint or to maximize a welfare/utility function. This a resource/market
allocation problem. PE/GE maximize consumer surplus, GR maximizes utility in a
Ramsey-Solow growth model.

PARTIAL EQUILIBRIUM (ENERGY SYSTEM), GENERAL EQUILIBRIUM OR OPTIMAL GROWTH:
---------------------------------------------------------------------------

    * MARKAL, TIMES PE (ETSAP)
    * MARKAL-MACRO, TIMES-MACRO GE (ETSAP)
    * EPPA recursive dynamic GE (MIT)
    * MESSAGE PE (fortran matrix generator + CLPEX??) (IIASA)
    * WITCH GR (FEEM)
    * REMIND-R GR (PIK)
    * DICE, RICE GR  (Nordhaus)
    * CRED GR (SEI)
    * MERGE GR (EPRI)


SIMULATION OR (discrete) SYSTEM DYNAMICS MODELS:
=====================================

In simulation or system dynamics model the allocation problem is solved by
'algorithm'.

ENERGY SYSTEM (bottom up or hybrid):
------------------------------------

    * NEMS(EIA)
    * WEPS+(EIA)
    * TIMER (PBL, RIVM)


INTEGRATED ASSESSMENT MODELS:
-----------------------------

    * GCAM - socioeconomic section (PNNL)




INTEGRATED ASSESSMENT MODELS WITH INTERMEDIATE COMPLEXITY CLIMATE MODELS:
=========================================================================

Note: The climate and socio-economic models are very loosely coupled. The
socio-economic part (macro+energy+land use) is run first. Emissions etc generated
are used to run the intermediate complexity climate model (MAGICC/ MIT models) and
climate change effects are calculated.