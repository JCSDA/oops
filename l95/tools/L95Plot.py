# (C) Copyright 2009-2016 ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

class L95Plot:

  #===========================================================================
  # __init__
  # Constructor:
  # Create the Figure for the plot
  #===========================================================================
  def __init__(self):
    self.__f = Figure()
    self.state = None
    self.increment = None  
    
  #=============================================================================
  # setPlot
  # put the datas to plot
  #=============================================================================
  def setPlot(self, title, expeData, bgData, obsData):
    n = len(expeData)
    x = np.arange(0, n)
    xobs = []

    self.title = self.__f.suptitle(title)
    self.state = self.__f.add_subplot(211)
    #self.state.set_title("State:")
    self.state.plot(x, expeData, 'b', label='state')
    
    if obsData is not None:
      for iobs in obsData['x']:
        xobs.append(n * iobs)
      yobs = obsData['y']
      self.state.plot(xobs, yobs, 'rx', label='obs')
      
    mini = np.amin(expeData) - 1
    maxi = np.amax(expeData) + 3
    self.state.axis([0, n, mini, maxi])
    
    # Compute increment Data
    incrData = None
    if bgData is not None:
      incrData = []
      for i in range(0, len(expeData)):
        incrData.append(expeData[i] - bgData[i])
      self.state.plot(x, bgData, 'b:', label='background')
    self.state.legend()

    if xobs is not None:
        yDep = []
        for i in range(0, len(xobs)):
          yDep.append(float(yobs[i]) - bgData[int(xobs[i])])      
              
    self.increment = self.__f.add_subplot(212)
    #self.increment.set_title("Increment:")
    self.increment.plot([0, n], [0.0, 0.0], 'k-', label='')

    ax = self.__f.gca()
    ax.set_xticks(np.arange(0, 40, 5))
    ax.set_yticks(np.arange(-5., 5., 0.5))
    self.increment.grid()
    if incrData is not None:
      self.increment.plot(x, incrData, 'b', label='increment')
    if len(xobs) > 0:
      self.increment.plot(xobs, yDep, 'rx', label='departure')
    self.increment.axis([0, n, -2, 2])
    self.increment.legend()

  #=============================================================================
  # getFigure
  # Accessor for the figure
  #=============================================================================
  def getFigure(self):
      return self.__f
  
  #=============================================================================
  # clear
  # Clear the plot
  #=============================================================================
  def clear(self):
      self.__f.texts = []
      if self.state is not None:
         self.state.clear()
      if self.increment is not None:
         self.increment.clear()
      self.state = None
      self.increment = None         
