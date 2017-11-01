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
import glob
import os
import re
import os.path
import xml.dom.minidom
from xml.dom.minidom import Node

class Experiment:

  #=============================================================================
  # Get some informations from the experiment's files found in the directory
  #=============================================================================
  def __init__(self, dir):
    self.dir = dir
    self.__map = {}
    l = glob.glob(dir + '/*.xml') 
    for i in l:
      filenameSplitted = os.path.split(i)[1].split('.')
      id = filenameSplitted[0]
      file = self.getFromXML(i)
      self.__map[id] = file

  #=============================================================================
  # getFromXML
  #=============================================================================      
  def getFromXML(self, ficXML):
    ret = {}
    dom = xml.dom.minidom.parse(ficXML)
    
    # Get output filename
    outputs = dom.getElementsByTagName('output')
    for output in outputs:
      ret['datadir'] = output.getElementsByTagName('datadir')[0].firstChild.data
      ret['exp'] = output.getElementsByTagName('exp')[0].firstChild.data
      ret['type'] = output.getElementsByTagName('type')[0].firstChild.data

    # Get background filenames
    bg = {}
    backgrounds = dom.getElementsByTagName('background')
    for background in backgrounds:
      states = background.getElementsByTagName('state')
      for state in states:
        bg['expe'], bg['refTime'] = self.getStateFromXML(state)
    initials = dom.getElementsByTagName('initial')
    for initial in initials:
      bg['expe'], bg['refTime'] = self.getStateFromXML(initial)
    ret['background'] = bg
    
    # Get obs filename
    obsFiles = dom.getElementsByTagName('ObsFileIn')
    ret['obsFileIn'] = None
    for obsFile in obsFiles:
      ret['obsFileIn'] = obsFile.firstChild.data
    return ret              
      
  #=============================================================================
  # getStateFromXML
  #=============================================================================         
  def getStateFromXML(self, node):
    bg = None
    filename = None
    try:
      filename = node.getElementsByTagName('filename')[0].firstChild.data
      date = node.getElementsByTagName('date')[0].firstChild.data
      f1 = filename.split('/')
      f2 = f1[len(f1) - 1]
      bg = f2.split('.')
      refTime = parser.parse(bg[2])      
    except: 
      return "xxx", "yyy"
    return bg[0] + "." + bg[1], refTime

      
  #=============================================================================
  # getIds
  # Returns an array containing all experiment's Ids
  #=============================================================================
  def getIds(self):
    return self.__map.keys()

  #=============================================================================
  # getObsFile
  #=============================================================================
  def getObsFile(self, id):
    return self.__map[id]['obsFileIn']

  #=============================================================================
  # getOutputFile
  # Returns filename without reference time
  #=============================================================================
  def getOutputFile(self, id):
    ret = None
    try:
      f = self.__map[id]
      if "datadir" in f:
        return f["exp"] + "." + f["type"]
    except ValueError, Argument:
      print Argument
    return ret

  #=============================================================================
  # getBackground
  # Returns the background: an associative array with expe.type and refTime as keys
  #=============================================================================
  def getBackground(self, id):
    return self.__map[id]['background']

