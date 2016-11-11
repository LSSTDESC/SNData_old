from distutils.core import setup
import sys
import os
import re

packageName = 'sndata'
packageDir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          packageName)
versionFile = os.path.join(packageDir, 'version.py')

print(versionFile, packageName, packageDir)
with open(versionFile, 'r') as f:
      s = f.read()
# Look up the string value assigned to __version__ in version.py using regexp
versionRegExp = re.compile("__version__ = \"(.*?)\"")
# Assign to __version__
__version__ =  versionRegExp.findall(s)[0]
setup(# package information
      name=packageName,
      version=__version__,
      description='A set of utilities to analyze SN light curves',
      long_description=''' ''',
      # What code to include as packages
      packages=[packageName],
      package_dir={packageName:'analyzeSN'},
      # What data to include as packages
      include_package_data=True,
      package_data={packageName: ['example_data/*.FITS', 'example_data/*.dat',
                                  'example_data/*.md']}
      )
