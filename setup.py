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
Installer for recem
"""

import os
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(name='Recem',
      version='0.1',
      maintainer='Shoibal Chakravarty',
      maintainer_email='shoibalc@gmail.com',
      url = 'http://shoibalc.github.com/Recem/',
      license = 'BSD',
      platforms = ["any"],
      packages = ['recem'],
      description = "Repository of Energy, Climate, and Economics Models",
      long_description = read('README.txt'),
      classifiers = [
            'Development Status :: 1 - Planning',
            'Intended Audience :: End Users/Desktop',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Natural Language :: English',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: Unix',
            'Programming Language :: Python',
            'Programming Language :: Unix Shell',
            'Topic :: Scientific/Engineering'
        ],
      )