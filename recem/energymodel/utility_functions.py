#==============================================================================
# RECEM: Repository of Energy, Climate, and Economics Models
#
#FIXME Copyright (c) 2012 Shoibal Chakravarty shoibalc[at]gmail[dot]com
#  [should it be Princeton University]
# This software is distributed under the BSD License.
# You may obtain a copy of the License at
# http://www.opensource.org/licenses/bsd-license.php
#==============================================================================



from coopr.pyomo import *


def copy_PAST(model):
    return model.PAST


def copy_HORIZON(model):
    return model.HORIZON


def copy_TIME(model):
    return model.TIME


def calculate_PERIOD(model, time):
    pos = model.HORIZON.ord(time)
    if pos == 1:
        return 0
    else:
        return (time - model.HORIZON[pos-1])


def calculate_PERIODDISCOUNT(model, time):

    init = model.HORIZON.first()
    years = time - init
    return ((1 - value(model.DISCOUNTRATE)) ** years)


def build_basic_model_utility(model):

# SETS and PARAMETERS [that have to do with time]

    model.PAST              = Set(ordered=True, within=Integers)
    model.HORIZON           = Set(ordered=True, within=Integers)
#        model.DATATIME     = Set(ordered=True, within=Integers)
    model.TIME              = model.PAST | model.HORIZON
    model.SEASON            = Set()
    model.TIMESLICE         = Set()
    model.TIMESLICESHARE    = Param(model.SEASON, model.TIMESLICE)
    model.PERIOD            = Param(model.HORIZON, rule=calculate_PERIOD)

# NOTE: In Python, B = A for complex objects [like model.PAST] only makes B
# refer to the A but does not make a copy of A. If we want a copy
# we have to create a new object B and copy the data of A. That is what we do
# in the next few lines. Functions copy_PAST, copy_FUTURE etc are defined in
# this module [recem.energymodel.utility_functions]

    model.VINTAGEPAST    = Set(ordered=True, initialize=copy_PAST)
    model.VINTAGEHORIZON = Set(ordered=True, initialize=copy_HORIZON)
    model.VINTAGETIME    = Set(ordered=True, initialize=copy_TIME)

# SETS [that have to do with demands, technology processes to meet these
#       demands, and the commodites that flow in and (flow) out of these
#       processes. We also list commodity imports and exports here. CO2 is
#       also a 'commodity' produced by technology processes!]

    model.DEMANDTYPES   = Set()
    model.TECH          = Set()
    model.BASELOAD      = Set(within=model.TECH)
    model.PEAKING       = Set(within=model.TECH)
    model.COMMODITY     = Set()
    model.GHG           = Set(within=model.COMMODITY, initialize=['CO2'])
    model.FLOW_IN       = Set(model.TECH,within=model.COMMODITY)
    model.FLOW_OUT      = Set(model.TECH,within=model.COMMODITY)
#    model.FLOW          = model.FLOW_IN.union(model.FLOW_OUT)
    model.IMPORT        = Set(model.COMMODITY)
    model.EXPORT        = Set(model.COMMODITY)

# PARAMETERS [discount rate: annual and for the length of individual periods]

    model.DISCOUNTRATE   = Param()
    model.PERIODDISCOUNT = Param(model.HORIZON, rule=calculate_PERIODDISCOUNT)

# PARAMETERS [demand, technoogy and flow related]

#    model.CURRENTCAPACITY   = Param(model.TECH, model.VINTAGEPAST)
#    model.EFFICIENCY        = \
#        Param(model.FLOW_IN, model.TECH, model.VINTAGEPAST, model.FLOW_OUT)