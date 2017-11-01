# (C) Copyright 2009-2016 ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

#===============================================================================
# QG0db
#
# A class to handle the observations used in QG model
#===============================================================================
import numpy as np
from dateutil import parser
from datetime import tzinfo, timedelta, datetime
import glob
import os
import re
import os.path

class QGFdb:

  #=============================================================================
  # Construct the fdb object, read the filenames and store them in a map
  #=============================================================================
  def __init__(self, fdbdir):
    self.fdbdir = fdbdir
    self.map = {}
    l = glob.glob(fdbdir + '/*.*.*') 
    for i in l:
      file = {}
      filenameSplitted = os.path.split(i)[1].split('.')
      if len(filenameSplitted[2]) < 10:
        continue
      if i.find(".obt") != -1:
        continue
      file["expe"] = filenameSplitted[0]
      file["type"] = filenameSplitted[1]
      file["referenceTime"] = filenameSplitted[2]
      file["refTime"] = parser.parse(filenameSplitted[2]) 
      try:
        term = self.getTerm(filenameSplitted[3])
        file["term"] = filenameSplitted[3]
        file["intTerm"] = term
        file["validTime"] = file["refTime"] + term

      except:
        file["intTerm"] = timedelta(seconds=0)
        file["validTime"] = file["refTime"]
        file["term"] = None
      file["filename"] = os.path.split(i)[1]
      self.map[os.path.split(i)[1]] = file
      
  #=============================================================================
  # getAtValidTime
  #=============================================================================
  def getAtValidTime(self, date):
      ret = []
      for file in self.map.values():
          if file["validTime"] == date :
              ret.append(file)
      return ret

  #=============================================================================
  # ggetAtValidTime
  #=============================================================================
  def ggetAtValidTime(self, stringdate):
      date = parser.parse(stringdate)
      return self.getAtValidTime(date)
  
  #=============================================================================
  # getData
  #=============================================================================
  def getData(self, file):
    if "psi" in file:
       psi = file["psi"]
       pv  = file["pv"] 
    else:
       infile = open(self.fdbdir + "/" + file["filename"], "r")
       
       [nx, ny, nl, nf, ns] = [int(i) for i in infile.readline().split()]
       validity_date = infile.readline()
       psi = np.fromfile(file=infile,dtype="float64", \
                      count=(2*ny*nx),sep=" ").reshape((2,ny,nx))
       pv  = np.fromfile(file=infile,dtype="float64", \
                      count=(2*ny*nx),sep=" ").reshape((2,ny,nx))
       infile.close()
                    
       file["psi"] = psi
       file["pv"] = pv
       
    ret = {}
    ret['pv'] = pv
    ret['psi'] = psi
    return ret

  #=============================================================================
  # getMetaData
  #=============================================================================
  def getMetaData(self):
    lst = {}
    for f in self.map.values():
        key = f["expe"] + "." + f["type"]
        if key in lst:
            arr = lst[key]
        else:
            arr = {}
            lst[key] = arr
        if f["referenceTime"] in arr:
            terms = arr[f["referenceTime"]]
        else:
            terms = {}
            arr[f["referenceTime"]] = terms
        terms[f["intTerm"]] = f["term"]
    for expe in lst.keys():
        for date in lst[expe].keys():
           t = lst[expe][date];
           k = t.keys()
           k.sort()
           newT = []
           for i in k:
               newT.append(t[i])
           lst[expe][date] = newT

    return lst

  #=============================================================================
  # getTimeList
  #=============================================================================
  def getTimeList(self, expe):
    lst = []
    items = {}
    for f in self.map.values():
      key = f["expe"] + "." + f["type"]
      if key == expe:
        items[f['refTime']] = f['referenceTime']
    k = items.keys()
    k.sort()
    for i in k:
      lst.append(items[i]) 
    return lst

  #=============================================================================
  # getExpeAt
  #=============================================================================
  def getExpeAt(self, expe, refTime, validTime):
    file = None
    for f in self.map.values():
      key = f["expe"] + "." + f["type"]
      if expe == key and refTime == f['refTime'] and validTime == f['validTime'] :
        file = f
    return file

  #=============================================================================
  # printMetaData
  #=============================================================================
  def printMetaData(self):
    lst = self.getMetaData()
    for expe in lst.keys():
        print expe
        for date in lst[expe].keys():
            print "-- ", date
            if len(lst[expe][date]) > 1:
              for ech in lst[expe][date]:
                print "---- ", ech
            

    
  #=============================================================================
  # getTerm
  #=============================================================================
  def getTerm(self, term):
    regex = re.compile('(?P<sign>-?)P(?:(?P<years>\d+)Y)?(?:(?P<months>\d+)M)?(?:(?P<days>\d+)D)?(?:T(?:(?P<hours>\d+)H)?(?:(?P<minutes>\d+)M)?(?:(?P<seconds>\d+)S)?)?')
    duration = regex.match(term).groupdict(0)
    delta = timedelta(days=int(duration['days']) + (int(duration['months']) * 30) + (int(duration['years']) * 365),
                  hours=int(duration['hours']),
                  minutes=int(duration['minutes']),
                  seconds=int(duration['seconds']))
    return delta
    
