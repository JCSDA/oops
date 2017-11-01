# (C) Copyright 2009-2016 ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

#===============================================================================
# L950db
#
# A class to handle the observations used in L95 model
#===============================================================================
import numpy as np
from dateutil import parser
from datetime import tzinfo, timedelta, datetime

class L95Odb:

  #=============================================================================
  # Construct the odb object, read the observations and store them in a map
  #=============================================================================
  def __init__(self, fdbdir, odb_file):
      self.datas = []
      infile = open(fdbdir + "/" + odb_file, "r")
      nbCols = int(infile.readline())
      columns = []
      for i in range(0, nbCols):
          columns.append(infile.readline().strip())
      nbRows = int(infile.readline())
      for i in range(0, nbRows):
          row = {}
          col = infile.readline().split();
          row["stringdate"] = col[1]
          d = parser.parse(col[1])
          row["date"] = d - timedelta(minutes=1)
          row["x"] = col[2]
          for j in range(0, nbCols) :
              row[columns[j]] = col[j + 3]
          self.datas.append(row)
          
  #=============================================================================
  # get the observations at a date
  #=============================================================================
  def get(self, date):
      obs = []
      for o in self.datas:
          if o["date"] == date:
              obs.append(o)
      return obs
  
  #=============================================================================
  # get the observations at a date
  #=============================================================================
  def getXYAt(self, date, colName):
      x = []
      y = []
      #try:
      for o in self.datas:
          if o["date"] == date:
            x.append(float(o["x"]))
            y.append(o[colName])
      #except: 
      #    x = []
      #    y = []
      return x,y
      
  #=============================================================================
  # get the observations at a date
  #=============================================================================
  def getAt(self, stringdate):
      date = parser.parse(stringdate)
      return self.get(date)
      