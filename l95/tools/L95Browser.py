# (C) Copyright 2009-2016 ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

from L95Plot import L95Plot
from L95Odb import L95Odb
from L95Fdb import L95Fdb

#===============================================================================
# L95Factory
# The factory for Lorenz 95 model
#===============================================================================
class L95Factory:

  #=============================================================================
  # createPlot
  #=============================================================================
  def createPlot(self):
    return L95Plot()

  #=============================================================================
  # createFdb
  #=============================================================================
  def createFdb(self, dir):
    return L95Fdb(dir)

  #=============================================================================
  # createOdb
  #=============================================================================
  def createOdb(self, file):
    return L95Odb(file)

  #=============================================================================
  # title
  #=============================================================================
  def title(self):
    return "Lorenz 95"

  #=============================================================================
  # id
  #=============================================================================
  def id(self):
    return "L95"

'''
  Program main
'''
if __name__ == "__main__":
    import sys
    import os

    workDir = os.path.split(__file__)[0]
    if len(workDir) > 0:
      os.chdir(workDir)
    sys.path.append("../../browser")

    from optparse import OptionParser
    from Browser import Browser

    parser = OptionParser()
    parser.add_option("-d", "--datadir", dest="dataDir", help="Data directory")
    parser.add_option("-c", "--configdir", dest="configDir", help="XML files directory")
    (options, args) = parser.parse_args()
    if options.configDir is None:
       options.configDir = "../config"
    if options.dataDir is None:
       options.dataDir = "../data"
    Browser(L95Factory(), options.configDir, options.dataDir)
