#==============================================================================
# RECEM: Repository of Energy, Climate, and Economics Models
#
# Copyright (c) 2012 Princeton University
# Author: Shoibal Chakravarty
# This software is distributed under the BSD License.
# You may obtain a copy of the License at
# http://www.opensource.org/licenses/bsd-license.php
#==============================================================================

__all__ = ['DICE_results_writer']
import sys
import time
import types
from coopr.pyomo import *

def DICE_results_writer(results, instance, file=None, mode='w', startyear=None):
    """
    DICE_results_writer is  functions that writes the results of solve process for the
    DICE models.
    NOTE: This is guaranteed to work for DICE models only!
    NOTE: Some columns report None which is equivalent to no value (or NaN).
    print_results(instance,file=FILE, mode=MODE) takes the relevant DICE
    instance, results returned by the solver and prints the
    1.  Solver summary
    2. Objective,
    3. Variables [lower bound, value, upper bound]
    4. Constraints [lower bound, value, upper bound, dual (if available)]
    in space delimited format. FILE is the name of the output file.
    MODE is write mode: 'w' (write) or 'a' (append). The defaults are:
    (FILE: results.txt, MODE: 'w')

    """
# NOTE: Except for the fact that we print the period in years instead of as an
# index, this function is very general. Should write up a standard version.
    instance.load(results)
    if file is not None:
        fp = open(file, mode)
    else:
        fp = open('results.txt', mode)
    if startyear is None:
        print "startyear is 1965 for DICE1993 and 2005 for DICE2007"
        return

    print >>fp, '\"',instance.name,'\"'
    print >>fp, '\"',instance.doc,'\"'
    print >>fp, '\"Solver Summary\"'
    print >>fp, '\"', results['Solver'][0], '\"'
    # Objective
    obj = instance.active_components(Objective)
    if len(obj) > 1:
        print >>fp, """
        Warning: More than one objective.  Using first objective.
        """
    print >>fp, "Objective\n", '\"', obj[obj.keys()[0]].doc, '\"'
    print >>fp, obj[obj.keys()[0]][None].expr()
    # Variables
    for v in instance.active_components(Var):
        varobject = getattr(instance, v)
        print >> fp, ""
        print >> fp, "Variable",v, '\"',varobject.doc,'\"'
        print >> fp, 'PERIOD  ',v, '  LOWER', 'VALUE', 'UPPER'
    # Special condition for singleton Variable (TRANS) which has no index
        if type(varobject.index) is types.NoneType:
            print>>fp, "-", v, varobject.lb is None and '-INF' or varobject.lb, \
            varobject.value, varobject.ub is None and '+INF' or varobject.ub
        else:
            for index in varobject:
                print >> fp, str(startyear+10*(index-1))+"-"+str(startyear+10*index),varobject[index].name,\
                varobject[index].lb is None and '-INF' or varobject[index].lb, \
                varobject[index].value, varobject[index].ub is None and '+INF' or varobject[index].ub

    # Constraints (duals, if available)
    print >> fp, "Constraints\n\n"
    for c in instance.active_components(Constraint):
        cobject = getattr(instance, c)
        print >> fp, ""
        print >> fp, "Constraint", c, '\"',cobject.doc,'\"'
        print >> fp, 'PERIOD  ', c, 'LOWER', 'VALUE', 'UPPER', 'DUAL'
    # Special condition for singleton Constraint (TRANSE) which has no index
        if type(cobject.index()) is dict and cobject.index().keys()[0] is None:
            print >> fp, "-", cobject.name,\
            cobject[None].lower is None and '-INF' or cobject[None].lower(),\
            cobject[None].body(), \
            cobject[None].upper is None and '+INF' or cobject[None].upper(),\
            cobject[None].dual
        else:
            for index in cobject:
                print >> fp, str(startyear+10*(index-1))+"-"+str(startyear+10*index), cobject[index].name,\
                cobject[index].lower is None and '-INF' or cobject[index].lower(),\
                cobject[index].body(), \
                cobject[index].upper is None and '+INF' or cobject[index].upper(),\
                cobject[index].dual