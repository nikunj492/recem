#==============================================================================
# RECEM: Repository of Energy, Climate, and Economics Models
#
#FIXME Copyright (c) 2012 Shoibal Chakravarty shoibalc[at]gmail[dot]com
#  [should it be Princeton University]
# This software is distributed under the BSD License.
# You may obtain a copy of the License at
# http://www.opensource.org/licenses/bsd-license.php
#==============================================================================

"""
EnergyModel is a bottom-up model in the TIMES/MARKAL
family of models. The EnergyModel is written using the AbstractModel
class of the Pyomo modeling language of the Coopr optimization tools
packages. See

https://software.sandia.gov/trac/coopr

The model is an abstract optimization model which is filled with data from
an appropriate data file (a *.dat file in the AMPL data format). This is then
translated by Pyomo to an instance file (another special data format file) and
sent to specialized solver programs like GLPK, CBC, IPOPT. The solution is
read by Pyomo and reported in human readable format.
"""

__all__ = ['EnergyModel']


from utility_functions import *
import logging
logger = logging.getLogger('recem')


class EnergyModel(AbstractModel):
    """
    The EnergyModel is an abstract energy-economics model to be filled
    with user supplied data. It is pre-packaged with a standard structure
    comes from the standard design and structure of most such models. The
    EnergyModel can be easily extended using the addBlock, addPolicy and
    similar methods.

    Technically, the EnergyModel class is a light wrapper  around the
    AbstractModel class of the Pyomo modeling language of the Coopr
    optimization tools packages [see coopr.pyomo.base.PyomoModel.AbstractModel]

    See, for further intformation on  Coopr:

    https://software.sandia.gov/trac/coopr

    "AbstractModel is an abstract optimization model that defers
    construction of components."
    """

    def __init__(self, *args, **kwds):
        AbstractModel.__init__(self, *args, **kwds)
        self.build_basic_model()
    pyutilib.component.core.alias("EnergyModel", 'EnergyModel is an\
    abstract energy optimization model')

    def build_basic_model(self):
        """
        Fill the AbstractModel with a MARKAL type abstract energymodel.
        method build_masic_model calls build_basic_model_utility in
        module utility_functions (see recem.energymodel.utility_functions)
        """
        build_basic_model_utility(self)