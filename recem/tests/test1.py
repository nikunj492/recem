#==============================================================================
# RECEM: Repository of Energy, Climate, and Economics Models
#
#FIXME Copyright (c) 2012 Shoibal Chakravarty shoibalc[at]gmail[dot]com
#  [should it be Princeton University]
# This software is distributed under the BSD License.
# You may obtain a copy of the License at
# http://www.opensource.org/licenses/bsd-license.php
#==============================================================================

from recem.energymodel import EnergyModel
from coopr.pyomo import *
import shutil

TestModel = EnergyModel("Testing!")
TestModel.build_basic_model()
shutil.copy("TestData.txt", "TestData.dat")
instance = TestModel.create('TestData.dat')
instance.pprint()