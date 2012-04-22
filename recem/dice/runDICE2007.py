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

from coopr.pyomo import *
from coopr.opt.base import SolverFactory
from coopr.opt.parallel import SolverManagerFactory
from recem.dice import createDICE2007
from recem.dice import DICE_results_writer
#from DICE2007 import createDICE2007
#from DICEutils import DICE_results_writer

global start_time
start_time = time.time()

dice = createDICE2007()
dice.doc = 'OPTIMAL SCENARIO'
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
print '[%8.2f] create model %s OPTIMAL SCENARIO\n' %(time.time()-start_time,dice.name)

instance = dice.create()
print '[%8.2f] created instance\n' %(time.time()-start_time)

results = solver_manager.solve(instance, opt=opt, tee=tee, options=options, suffixes=['dual','rc'])
print '[%8.2f] Solved instance\n' %((time.time()-start_time))

DICE_results_writer(results, instance, file='dice2007_optimum.txt', mode='w', startyear=2005)
print ""
print '[%8.2f] Finished writing results.\n' %((time.time()-start_time))

