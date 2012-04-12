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

import os
import sys
import time

from coopr.pyomo import *
from coopr.opt.base import SolverFactory
from coopr.opt.parallel import SolverManagerFactory

global start_time
start_time = time.time()

def printvaluesbyyear(T, R, Name="Emissions (GtC/decade)"):
    print "period      ",Name
    for t in T:
        if t is not T.last():
            print str(1965+10*(t-1))+"-"+str(1965+10*t),"  ", value(R[t])


from recem.dice import createDICE1993
opt = SolverFactory('ipopt',solver_io='nl')
options =   """
        halt_on_ampl_error=yes
        constr_viol_tol=1e-1
        acceptable_compl_inf_tol=1e-3
        acceptable_dual_inf_tol=1e+10
        tol=1e-4
        print_level=4
        max_iter=4000
        """
tee = False

solver_manager = SolverManagerFactory('serial')
#modeldata = ModelData()
#modeldata.add('dice-v3.dat')

dice = createDICE1993()

# //////////===== BASE SCENARIO =====\\\\\\\\\\
# In the base scenario, there is no abatement and hence no abatement costs
# Moreover, the temperature costs are also assumed to be zero (implemented through BASE=0).

# In the base, there is no abatement, so MIU has to be equal to zero.

#fix {t in T} MIU[t] := 0;
dice.MIU.bounds = (0,0)   # fix by setting upper and lower bounds to 0.
#dice.MIU.fixed = True
#let BASE := 0;
dice.BASE.set_value(0)

print '[%8.2f] create model BASE SCENARIO\n' %(time.time()-start_time)

print """In the base scenario, there is no abatement
 and hence no abatement costs Moreover, the temperature costs
 are also assumed to be zero (implemented through BASE=0)."""

inst = dice.create()

print '[%8.2f] created instance\n' %(time.time()-start_time)

results = solver_manager.solve(inst,opt=opt,tee=tee, options = options)
print '[%8.2f] Solved instance\n' %((time.time()-start_time))

transformed_results = inst.update_results(results)
#transformed_results.write(filename='BASE.yml')
#print '[%8.2f] Write results to file.\n' %((time.time()-start_time))

inst.load(results)
print "Emissions: Base scenario"
#display(inst.MIU)
printvaluesbyyear(T=inst.T, R=inst.E, Name="Emissions (GtC/decade)")
#display(inst.OBJ)
#display(inst.KK)
print ""
#print '[%8.2f] Finished summary.\n' %((time.time()-start_time))

# //////////===== MARKET SCENARIO =====\\\\\\\\\\
# In this scenario, there are temperature costs (BASE = 1)
# but as abatement is still fixed at zero, there are no abatement costs.
# Nordhaus has labelled this scenario 'the market solution'.

#SCENARIO(VERSIONS) = NO;
#SCENARIO('MARKET') = YES;
#BASE        = 1;

#fix {t in T} MIU[t] := 0;
dice.MIU.bounds = (0,0)   # fix by setting upper and lower bounds to 0.

#let BASE := 1;
dice.BASE.set_value(1)

print '[%8.2f] create model MARKET SCENARIO\n' %(time.time()-start_time)
print """In this scenario, there are temperature costs (BASE = 1)
    but as abatement is still fixed at zero, there are no abatement costs.
    Nordhaus has labelled this scenario 'the market solution'."""

inst = dice.create()
#inst.write(filename='prob.nl',format='nl')

print '[%8.2f] created instance\n' %(time.time()-start_time)

results = solver_manager.solve(inst,opt=opt,tee=tee)
print '[%8.2f] Solved instance\n' %((time.time()-start_time))

transformed_results = inst.update_results(results)
#transformed_results.write(filename='MARKET.yml')
#print '[%8.2f] Write results to file.\n' %((time.time()-start_time))

inst.load(results)
print "Emissions: Market scenario"
#display(inst.MIU)
printvaluesbyyear(T=inst.T, R=inst.E, Name="Emissions (GtC/decade)")
#display(inst.OBJ)
#display(inst.KK)
print ""
#print '[%8.2f] Finished summary.\n' %((time.time()-start_time))




#//////////===== OPT_CONT SCENARIO =====\\\\\\\\\\
# In this third scenario, the optimal abatement rate is calculated.
# This optimal point is where marginal abatement costs equal marginal temperature costs.

#SCENARIO(VERSIONS) = NO;
#SCENARIO('OPT_CONT') = YES;

# The variable MIU is 'freed' by the .UP and .LO statements.
# Abatement in the (historical) first three periods is fixed at zero.

#MIU.UP(T)   = 0.99;     MIU.LO(T)   = 0.01;
#MIU.FX('1') = 0;        MIU.FX('2') = 0;        MIU.FX('3') = 0;

#unfix {t in T} MIU[t] := 0;
dice.MIU.bounds = (0.0,0.99)   #unfix by restoring original lower and upper bounds.
dice.MIU._initialize = 0.005  # this shouldn't be done but I don't have another option:)
#dice.MIU = Var(dice.T, bounds=(0.0,0.99), initialize=0.05)
#fix MIU[1] := 0;
#fix MIU[2] := 0;
#fix MIU[3] := 0;
def define(dice, t):
    if 1<= t <= 3:
        return (0.0000001, dice.MIU[t], 0.0000001)
    else:
        return Constraint.Skip
dice.FIXMIU = Constraint(dice.T, rule = define) #return a tuple with equal up, lo values to fix vars.
#let BASE := 1;
dice.BASE.set_value(1)



print '[%8.2f] create model OPT_CONT SCENARIO\n' %(time.time()-start_time)
print """In this third scenario, the optimal abatement rate is calculated.
This optimal point is where marginal abatement costs equal marginal temperature costs."""


inst = dice.create()

print '[%8.2f] created instance\n' %(time.time()-start_time)

results = solver_manager.solve(inst,opt=opt,tee=tee)
print '[%8.2f] Solved instance\n' %((time.time()-start_time))

#results.write(filename='cbcnl.yml')
transformed_results = inst.update_results(results)
#transformed_results.write(filename='OPT_CONT.yml')
#print '[%8.2f] Write results to file.\n' %((time.time()-start_time))

inst.load(results)
print "Emissions: OTP_CONT Scenario"
#display(inst.MIU)
printvaluesbyyear(T=inst.T, R=inst.E, Name="Emissions (GtC/decade)")
#display(inst.OBJ)
#display(inst.KK)
print ""
#print '[%8.2f] Finished summary.\n' %((time.time()-start_time))


#//////////===== CONCENT SCENARIO =====\\\\\\\\\\
# In the last scenario, an upper bound is placed on CO2 concentrations.
# Note that the abatement rate is still free.

#SCENARIO(VERSIONS) = NO;
#SCENARIO('CONCENT') = YES;
#M.UP(T) = 1180;

#unfix {t in T} MIU[t] := 0;
#dice.MIU.bounds = (0.0001,0.99)   #unfix by restoring original lower and upper bounds.
#dice.FIXMIU.deactivate()

#let BASE := 1;
dice.BASE.set_value(1)

#redeclare var	M{T}	>=0 <=1180:=M0	;	#	set upper bound on CO2EQ concentration (billion ton)
#
dice.M.bounds = (0.0, 1180)

print '[%8.2f] create model CONCENT SCENARIO\n' %(time.time()-start_time)
print """In the last scenario, an upper bound is placed on CO2 concentrations.
Note that the abatement rate is still free."""

inst = dice.create()
#modeldata.read(dice)
#inst = dice.create(modeldata)
#inst.write(filename='prob.nl',format='nl')
print '[%8.2f] created instance\n' %(time.time()-start_time)

#outfile = 'dice-pyomo.nl'
#inst.write(filename=outfile,format='nl')
#sys.stdout.write('[%8.2f] wrote instance to %s\n' %((time.time()-start_time), outfile)


results = solver_manager.solve(inst,opt=opt,tee=tee)
print '[%8.2f] Solved instance\n' %((time.time()-start_time))

#results.write(filename='cbcnl.yml')
transformed_results = inst.update_results(results)
#transformed_results.write(filename='CONCENT.yml')
#sys.stdout.write('[%8.2f] Write results to file.\n' %((time.time()-start_time)))

inst.load(results)
print "Emissions: CONCENT Scenario"
#display(inst.MIU)
printvaluesbyyear(T=inst.T, R=inst.E, Name="Emissions (GtC/decade)")
#display(inst.OBJ)
#display(inst.KK)
print ""
sys.stdout.write('[%8.2f] Finished.\n' %((time.time()-start_time)))


