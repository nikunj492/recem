#==============================================================================
# RECEM: Repository of Energy, Climate, and Economics Models
#
# Copyright (c) 2012 Princeton University
# Author: Shoibal Chakravarty
# This software is distributed under the BSD License.
# You may obtain a copy of the License at
# http://www.opensource.org/licenses/bsd-license.php
#==============================================================================

EnergyModel README [TIMES/MARKAL/MESSAGE family of models]

BASIC SETS:
-----------

    r           REGION
    t,v         TIME, VINTAGE
    s           TIMESLICE
    p           PROCESS or TECHNOLOGY
    c           COMMODITY [types: energy, material, emissions, demand]
    in/out      BINARY FLAG SET [storage: in OR out]
    imp/exp     BINARY FLAG SET [trade: imp OR exp]

ENERGY VARIABLES:
-----------------

    NewCapacity[r,v,p]:
    Capacity Additions of PROCESS p in REGION r at TIME v

    Capacity[r,v,t,p]:
    Capacity of VINTAGE v for PROCESS p in REGION r at TIME t

    CapacityTotal[r,t,p]:
    Total Capacity for PROCESS p in REGION r at TIME t (depends on Capacity[r,v,t,p])

    Activity[r,v,t,p,s]:
    Activity of PROCESS p of VINTAGE v in REGION r at TIMESLICE s at TIME t

    Flow[r,v,t,p,c,s]:
    Quantity of COMMODITY c in REGION r consumed/produced by PROCESS p of VINTAGE v
    at TIMESLICE s at TIME t

    Storage[r,v,t,p,c,s,in/out]:
    Quantity of COMMODITY c stored[in]/discharged[out] in REGION r by PROCESS p of VINTAGE v
    at TIMESLICE s at TIME t

    Trade[r,t,p,c,s,imp/exp]:
    Quantity of COMMODITY c imported[imp]/exported[exp] to/from REGION r by PROCESS p
    at TIMESLICE s at TIME t

    Demand[r,t,d]:
    Demand of COMMODITY [of type demand] d in REGION r at TIME t

-------------------------------------------------------------------------------
#TODO CLIMATE CHANGE VARIABLES:
-------------------------------

    Concentration[t]
    ForcingFactor[t]
    TempAtmopshere[t]
    TempOcean[t]
    SeaLevelRise[r,t]

#TODO MACROECONOMIC VARIABLES:
------------------------------

    Population[r,t]
    Labor[r,t]
    Capital[r,t]
    Consumption[r,t]
-------------------------------------------------------------------------------

OBJECTIVE:
----------
