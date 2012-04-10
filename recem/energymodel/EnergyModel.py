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
family of models. The EnergyModel class is a light wrapper
around the AbstractModel class of the Pyomo modeling language
of the Coopr optimization tools packages. See

https://software.sandia.gov/trac/coopr

"""

__all__ = ['EnergyModel']

from coopr.pyomo import *
import pyutilib.component.core

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
    optimization tools packages [coopr.pyomo.base.PyomoModel.AbstractModel].
    See

    https://software.sandia.gov/trac/coopr

    AbstractModel is an abstract optimization model that defers
    construction of components.
    """

    def __init__(self, *args, **kwds):
        AbstractModel.__init__(self, *args, **kwds)

    pyutilib.component.core.alias("EnergyModel", 'EnergyModel is an\
    abstract energy optimization model')

    def bau(self):
        """

        """
#

#    def report(self):
#        """
#        Writes a tabulated and formatted report which is
#        friendly to spreadsheet programs. The content is
#        similar to the results file (results.yml or results.json).
#        """
#        pass
#
#    def plot(self):
#        pass