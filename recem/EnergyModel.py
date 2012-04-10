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
The EnergyModel is a bottom-up model in the TIMES/MARKAL
family of models.
"""

__all__ = ['EnergyModel']

from coopr.pyomo import *
import pyutilib.component.core

import logging
logger = logging.getLogger('recem.sem')


class EnergyModel(AbstractModel):
    """
    An abstract optimization model that defers construction of
    components.
    """

    def __init__(self, *args, **kwds):
        AbstractModel.__init__(self, *args, **kwds)

    pyutilib.component.core.alias("EnergyModel", 'EnergyModel is an\
    AbstractModel: An abstract optimization model that\
    defers construction of components.')