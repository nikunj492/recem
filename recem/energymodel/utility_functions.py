#==============================================================================
# RECEM: Repository of Energy, Climate, and Economics Models
#
# Copyright (c) 2012 Princeton University
# Author: Shoibal Chakravarty
# This software is distributed under the BSD License.
# You may obtain a copy of the License at
# http://www.opensource.org/licenses/bsd-license.php
#==============================================================================



from coopr.pyomo import *


def copy_PAST(m):
    return m.PAST


def copy_HORIZON(m):
    return m.HORIZON


def copy_TIME(m):
    return m.TIME


def calculate_PERIOD(m, time):
    pos = m.HORIZON.ord(time)
    if pos == 1:
        return 0
    else:
        return (time - m.HORIZON[pos-1])


def calculate_PERIODDISCOUNT(m, time):

    init = m.HORIZON.first()
    years = time - init
    return ((1 - value(m.DISCOUNTRATE)) ** years)


def build_basic_model_utility(m):

# SETS and PARAMETERS [that have to do with time]

    m.PAST              = Set(ordered=True, within=Integers)
    m.HORIZON           = Set(ordered=True, within=Integers)
#        m.DATATIME     = Set(ordered=True, within=Integers)
    m.TIME              = m.PAST | m.HORIZON
    m.SEASON            = Set()
    m.TIMESLICE         = Set()
    m.TIMESLICESHARE    = Param(m.SEASON, m.TIMESLICE)
    m.PERIOD            = Param(m.HORIZON, rule=calculate_PERIOD)

# NOTE: In Python, B = A for complex objects [like m.PAST] only makes B
# refer to the A but does not make a copy of A. If we want a copy
# we have to create a new object B and copy the data of A. That is what we do
# in the next few lines. Functions copy_PAST, copy_FUTURE etc are defined in
# this module [recem.energym.utility_functions]

    m.VINTAGEPAST    = Set(ordered=True, initialize=copy_PAST)
    m.VINTAGEHORIZON = Set(ordered=True, initialize=copy_HORIZON)
    m.VINTAGETIME    = Set(ordered=True, initialize=copy_TIME)

# SETS [that have to do with demands, technology processes to meet these
#       demands, and the commodites that flow in and (flow) out of these
#       processes. We also list commodity imports and exports here. CO2 is
#       also a 'commodity' produced by technology processes!]

    m.DEMANDTYPES   = Set()
    m.TECH          = Set()
    m.BASELOAD      = Set(within=m.TECH)
    m.PEAKING       = Set(within=m.TECH)
    m.COMMODITY     = Set()
    m.GHG           = Set(within=m.COMMODITY, initialize=['CO2'])
    m.FLOW_IN       = Set(m.TECH,within=m.COMMODITY)
    m.FLOW_OUT      = Set(m.TECH,within=m.COMMODITY)
#    m.FLOW          = m.FLOW_IN.union(m.FLOW_OUT)
    m.IMPORT        = Set(m.COMMODITY)
    m.EXPORT        = Set(m.COMMODITY)

# PARAMETERS [discount rate: annual and for the length of individual periods]

    m.DISCOUNTRATE   = Param()
    m.PERIODDISCOUNT = Param(m.HORIZON, rule=calculate_PERIODDISCOUNT)

# PARAMETERS [demand, technoogy and flow related]

#    m.CURRENTCAPACITY   = Param(m.TECH, m.VINTAGEPAST)
#    m.EFFICIENCY        = \
#        Param(m.FLOW_IN, m.TECH, m.VINTAGEPAST, m.FLOW_OUT)