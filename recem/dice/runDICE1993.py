#==============================================================================
# RECEM: Repository of Energy, Climate, and Economics Models
#
#   NOTE: DICE and all related programs in this package [recem.dice] are in
#   the public domain. DICE and various versions of it have been in the public
#   domain since the first version of the models were published by
#   William Nordhaus in 1979.
#
#   You may use this program at your own risk!
#==============================================================================

import sys
import time
import types

from coopr.pyomo import *
from coopr.opt.base import SolverFactory
from coopr.opt.parallel import SolverManagerFactory
from recem.dice import createDICE1993
from recem.dice import DICE_results_writer
#from DICE1993 import createDICE1993
#from DICEutils import DICE_results_writer

global start_time
start_time = time.time()

opt = SolverFactory('ipopt',solver_io='nl')
options =   """
        halt_on_ampl_error=yes"""
#        constr_viol_tol=1e-1
#        acceptable_compl_inf_tol=1e-3
#        acceptable_dual_inf_tol=1e+10
#        tol=1e-4
#        print_level=4
#        max_iter=4000
#        """
tee = False

solver_manager = SolverManagerFactory('serial')

dice = createDICE1993()

#    * Note that there is no equation to determine the optimal value of MIU(T).
#    * This variable is determined endogenously by GAMS as a 'free variable'.
#
#    * In the base, there is no abatement, so MIU has to be equal to zero.
#    MIU.FX(T)   = 0;


#    * //////////===== BASE SCENARIO =====\\\\\\\\\\
#    * In the base scenario, there is no abatement and hence no abatement costs
#    * Moreover, the temperature costs are also assumed to be zero (implemented through BASE=0).
#
#    SCENARIO(VERSIONS) = NO;
#    SCENARIO('BASE') = YES;
#    BASE = 0;
#    MIU.FX(T) = 0;
#
#    SOLVE DICE MAXIMIZING UTILITY USING NLP;

dice.MIU.bounds = (0.0,0.0)   # fix by setting upper and lower bounds to 0.
dice.BASE.set_value(0.0)

print '[%8.2f] Create model BASE SCENARIO\n' %(time.time()-start_time)

dice.doc =  """
    BASE SCENARIO
    In the base scenario, there is no abatement
    and hence no abatement costs Moreover, the temperature costs
    are also assumed to be zero (implemented through BASE=0).
    """

print dice.doc
inst = dice.create()

print '[%8.2f] Created instance\n' %(time.time()-start_time)

results = solver_manager.solve(inst, opt=opt, tee=tee, options=options, suffixes=['dual','rc'])

print '[%8.2f] Solved instance\n' %((time.time()-start_time))

inst.load(results)
DICE_results_writer(results, inst, file='dice1993.txt', mode='w', startyear=1965)
print ""

#    *$ontext
#    * //////////===== MARKET SCENARIO =====\\\\\\\\\\
#    * In this scenario, there are temperature costs (BASE = 1)
#    * but as abatement is still fixed at zero, there are no abatement costs.
#    * Nordhaus has labelled this scenario 'the market solution'.
#
#    SCENARIO(VERSIONS) = NO;
#    SCENARIO('MARKET') = YES;
#    BASE        = 1;
#
#    SOLVE DICE MAXIMIZING UTILITY USING NLP;

dice.MIU.bounds = (0.0,0.0)   # fix by setting upper and lower bounds to 0.
dice.BASE.set_value(1.0)

print '[%8.2f] Create model MARKET SCENARIO\n' %(time.time()-start_time)
dice.doc = """
    MARKET SCENARIO
    In this scenario, there are temperature costs (BASE = 1)
    but as abatement is still fixed at zero, there are no abatement costs.
    Nordhaus has labelled this scenario 'the market solution'.
    """
print dice.doc
inst = dice.create()

print '[%8.2f] Created instance\n' %(time.time()-start_time)

results = solver_manager.solve(inst, opt=opt, tee=tee, options=options, suffixes=['dual','rc'])
print '[%8.2f] Solved instance\n' %((time.time()-start_time))

inst.load(results)
DICE_results_writer(results, inst, file='dice1993.txt', mode='a', startyear=1965)
print ""


#    *//////////===== OPT_CONT SCENARIO =====\\\\\\\\\\
#    * In this third scenario, the optimal abatement rate is calculated.
#    * This optimal point is where marginal abatement costs equal marginal temperature costs.
#
#
#    SCENARIO(VERSIONS) = NO;
#    SCENARIO('OPT_CONT') = YES;
#
#    * The variable MIU is 'freed' by the .UP and .LO statements.
#    * Abatement in the (historical) first three periods is fixed at zero.
#    MIU.UP(T)   = 0.99;     MIU.LO(T)   = 0.01;
#    MIU.FX('1') = 0;        MIU.FX('2') = 0;        MIU.FX('3') = 0;
#
#    SOLVE DICE MAXIMIZING UTILITY USING NLP;

dice.MIU.bounds = (0.0,0.9999)
def FIXMIU_rule(dice, t):
    if 1<= t <= 3:
        return (0.00000001, dice.MIU[t], 0.00000001) #return a tuple with equal up, lo values to fix vars.
    else:
        return Constraint.Skip
dice.FIXMIU = Constraint(dice.T, rule = FIXMIU_rule)
dice.BASE.set_value(1.0)

print '[%8.2f] Create model OPT_CONT SCENARIO\n' %(time.time()-start_time)
dice.doc= """
    OPT_CONT SCENARIO
    In this third scenario, the optimal abatement
    rate is calculated. This optimal point is where
    marginal abatement costs equal marginal temperature costs.
    """

print dice.doc
inst = dice.create()
print '[%8.2f] Created instance\n' %(time.time()-start_time)
results = solver_manager.solve(inst, opt=opt, tee=tee, options=options, suffixes=['dual','rc'])
print '[%8.2f] Solved instance\n' %((time.time()-start_time))

inst.load(results)
DICE_results_writer(results, inst, file='dice1993.txt', mode='a', startyear=1965)
print ""

#    *//////////===== CONCENT SCENARIO =====\\\\\\\\\\
#    * In the last scenario, an upper bound is placed on CO2 concentrations.
#    * Note that the abatement rate is still free.
#
#    SCENARIO(VERSIONS) = NO;
#    SCENARIO('CONCENT') = YES;
#    M.UP(T) = 1180;
#
#    SOLVE DICE MAXIMIZING UTILITY USING NLP;

dice.BASE.set_value(1.0)
dice.M.bounds = (0.0, 1180)

print '[%8.2f] Create model CONCENT SCENARIO\n' %(time.time()-start_time)
dice.doc =  """
    CONCENT SCENARIO
    In the last scenario, an upper bound is placed on CO2 concentrations.
    Note that the abatement rate is still free.
    """

print dice.doc
inst = dice.create()
print '[%8.2f] Created instance\n' %(time.time()-start_time)

results = solver_manager.solve(inst, opt=opt, tee=tee, options=options, suffixes=['dual','rc'])
print '[%8.2f] Solved instance\n' %((time.time()-start_time))

inst.load(results)
DICE_results_writer(results, inst, file='dice1993.txt', mode='a', startyear=1965)
print ""
print '[%8.2f] Finished.\n' %((time.time()-start_time))