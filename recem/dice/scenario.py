#import gc
#import logging
import os
import sys
#import textwrap
#import traceback
#import types
import time

from coopr.pyomo import *
from coopr.opt import ProblemFormat
from coopr.opt.base import SolverFactory
from coopr.opt.parallel import SolverManagerFactory
import pyutilib.misc
#from pyutilib.component.core import ExtensionPoint, Plugin, implements
#from pyutilib.services import TempfileManager

global start_time
start_time = time.time()

from DICE2007_short import M
opt = SolverFactory('ipopt',solver_io='nl')
tee = True
options =   """
        halt_on_ampl_error=yes
        """
#        constr_viol_tol=1e-5
#        acceptable_compl_inf_tol=1e-3
#        acceptable_dual_inf_tol=1e-3
#        tol=1e-4
#        max_iter=500
#        """
solver_manager = SolverManagerFactory('serial')
#modeldata = ModelData()

#
#//////////////////// OPTIMAL SCENARIO  /////////////////////////////////
#
sys.stdout.write('[%8.2f] create model OPTIMAL SCENARIO\n' %(time.time()-start_time))
#ep = ExtensionPoint(IPyomoScriptCreateModelData)
#if len(ep)==1:
#    modeldata = ep.service().apply(model=model)
#else:
#    modeldata = ModelData()
#modeldata.read(dice)
inst = M.create()
#inst.write(filename='prob.nl',format='nl')
sys.stdout.write('[%8.2f] created instance\n' %(time.time()-start_time))

results = solver_manager.solve(inst,opt=opt,tee=tee, options = options)
sys.stdout.write('[%8.2f] Solved instance\n' %((time.time()-start_time)))

#results.write(filename='cbcnl.yml')
transformed_results = inst.update_results(results)
transformed_results.write(filename='OPT.yml')
sys.stdout.write('[%8.2f] Write results to file.\n' %((time.time()-start_time)))

inst.load(results)
print ""
display(inst.MIU)
display(inst.E)
display(inst.UTIL)
print ""

sys.stdout.write('[%8.2f] Finished summary.\n' %((time.time()-start_time)))

