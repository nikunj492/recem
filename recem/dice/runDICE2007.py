#==============================================================================
# RECEM: Repository of Energy, Climate, and Economics Models
#
#   NOTE: DICE and all related programs in this package [recem.dice] are in
#   the public domain. DICE and various versions of it have been in the public
#   domain since the first version of the models were published by
#   William Nordhaus in 1979.
#
#   No license whatsover!! You use this program at your own risk!
#==============================================================================
import os
import sys
import time

from coopr.pyomo import *
from coopr.opt import ProblemFormat
from coopr.opt.base import SolverFactory
from coopr.opt.parallel import SolverManagerFactory
from recem.dice import createDICE2007
import pyutilib.misc

global start_time
start_time = time.time()


dice = createDICE2007()

opt = SolverFactory('ipopt',solver_io='nl')
tee = False
options =   """
        halt_on_ampl_error=yes"""
#        constr_viol_tol=1e-5
#        acceptable_compl_inf_tol=1e-3
#        acceptable_dual_inf_tol=1e-3
#        tol=1e-4
#        max_iter=2000
#        """
solver_manager = SolverManagerFactory('serial')

#
#//////////////////// OPTIMAL SCENARIO  /////////////////////////////////
#
sys.stdout.write('[%8.2f] create model OPTIMAL SCENARIO\n' %(time.time()-start_time))

instance = dice.create()
sys.stdout.write('[%8.2f] created instance\n' %(time.time()-start_time))

results = solver_manager.solve(instance,opt=opt,tee=tee, options = options)
sys.stdout.write('[%8.2f] Solved instance\n' %((time.time()-start_time)))

transformed_results = instance.update_results(results)
transformed_results.write(filename='OPT.yml')
sys.stdout.write('[%8.2f] Write results to file.\n' %((time.time()-start_time)))

instance.load(results)
print "Emissions: Optimal Scenario"
display(instance.MIU)
display(instance.E)
display(instance.UTIL)
print ""

sys.stdout.write('[%8.2f] Finished summary.\n' %((time.time()-start_time)))

