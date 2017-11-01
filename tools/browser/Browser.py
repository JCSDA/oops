#!/usr/bin/env python
# (C) Copyright 2009-2016 ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


import glob
import os
import os.path
import random
import sys
import time
from datetime import tzinfo, timedelta, datetime
from dateutil import parser
import gobject

try:
    import pygtk
    pygtk.require("2.0")
except:
    pass
try:
    import gtk
    import gtk.glade
except:
    sys.exit(1)

from matplotlib.backends.backend_gtkcairo import FigureCanvasGTKCairo as FigureCanvas
from matplotlib.figure import Figure

from ModelFactory import ModelFactory
from Experiment import Experiment

class Browser:

    #===========================================================================
    # __init__
    # Constructor:
    # Load XML which describes the widgets tree, connects signals and
    # initializes combo-boxes
    #===========================================================================
    def __init__(self, factory, expDir, fdbDir):
        # Initialize 
        self.readyForUpdate = False
        self.plt = None
        self.background = None
        self.fdbdir = fdbDir
        self.expdir = expDir
        self.model = factory.id()
        ModelFactory.addFactory(factory)

        # Build the GUI
        gladefile = __file__.replace(".py", ".xml")
        self.wTree = gtk.glade.XML(gladefile)
        dic = { "on_quit_clicked" : self.onQuitClicked,
                "on_png_clicked" : self.onPngClicked,
                "onSelectExpe": self.onSelectExpe,
                "onSelectDate": self.onSelectDate,
                "onSelectTerm": self.onSelectTerm,
                "onToggleAnim": self.onToggleAnim,
            }
        self.wTree.signal_autoconnect(dic)
        window = self.wTree.get_widget("window")
        window.set_title(ModelFactory.title(self.model) + " experiment's browser")
        window.connect("destroy", self.onQuitClicked)
        self.animToggle = self.wTree.get_widget("toggleAnim")

        # Initialize the experience's list
        self.nDates = 1
        self.nTerms = 1

        # Fill in experiment's menu
        comboexpe = self.wTree.get_widget("comboexpe")
        self.experiments = Experiment(self.expdir)
        self.fdb = ModelFactory.createFdb(self.model, self.fdbdir)
        comboexpe.remove_text(0)
        for expe in self.experiments.getIds():
          comboexpe.append_text(expe)

        # Build a canvas for plotting
        self.plt = ModelFactory.createPlot(self.model)
        self.canvas = FigureCanvas(self.plt.getFigure())
        self.box = self.wTree.get_widget("eventbox1")
        self.box.add(self.canvas)
        self.canvas.show()

        # Open the GUI
        self.readyForUpdate = True
        comboexpe.set_active(0)
        gtk.main()


    #===========================================================================
    # onQuitClicked:
    # Quit the application
    #===========================================================================
    def onQuitClicked(self, widget):
        sys.exit()

    #===========================================================================
    # onPngClicked:
    # Save current plot in a PNG file 
    #===========================================================================
    def onPngClicked(self, widget):
        try:
          now = time.mktime(time.localtime())
          strNow = time.strftime("%Y-%m-%d-%H:%M:%S", time.gmtime(now))
          pngFile = os.environ["HOME"] + '/Images/' + strNow + '.png'
          self.plt.getFigure().savefig(pngFile)
        except ValueError, Argument:
          print "Error in export pdf:", Argument


    #===========================================================================
    # onToggleAnim:
    #===========================================================================
    def anim(self):
        radio = self.wTree.get_widget("animRefTime")
        if radio.get_active() is True:
          combo = self.wTree.get_widget("combodate")
          n = self.nDates
        else:
          combo = self.wTree.get_widget("comboterm")
          n = self.nTerms
          
        if n > 1:
          num = combo.get_active()
          combo.set_active(num + 1)
        return self.animToggle.get_active()

    def onToggleAnim(self, widget):
        if widget.get_active() is True:
          self.timer = gobject.timeout_add(1000, self.anim)

    #===========================================================================
    # onSelectExpe:
    # Select experiment in menu
    #===========================================================================
    def onSelectExpe(self, widget):
       if not self.readyForUpdate:
         return
       try:
          # Remove the date's and term's lists
          self.readyForUpdate = False
          combodate = self.wTree.get_widget("combodate")
          oldDate = combodate.get_active_text()
          for i in range(0, self.nDates):
            combodate.remove_text(0)
          self.nDates = 0

          # Read observations's file
          self.idExpe = widget.get_active_text()
          if self.idExpe is None:
            self.abort()
            return

          self.odb = None
          obsFile = self.experiments.getObsFile(self.idExpe)

          if obsFile is not None:
            self.odb = ModelFactory.createOdb(self.model, self.fdbdir + '/' + os.path.split(obsFile)[1])
          self.output = self.experiments.getOutputFile(self.idExpe)

          # Fill in date's menu
          numOldDate = 0
          for d in self.fdb.getTimeList(self.output):
            combodate.append_text(d)
            if d == oldDate:
               numOldDate = self.nDates
            self.nDates += 1
          self.readyForUpdate = True
          combodate.set_active(numOldDate)
       except ValueError, Argument:
         self.abort()
         print "Error fired here !!", Argument

    #===========================================================================
    # onSelectDate:
    # Select date in menu
    #===========================================================================
    def onSelectDate(self, widget):
       if not self.readyForUpdate:
         return
       try:
          # Remove term's menu
          self.readyForUpdate = False
          comboexpe = self.wTree.get_widget("comboexpe")
          comboterm = self.wTree.get_widget("comboterm")
          oldTerm = comboterm.get_active_text()

          for i in range(0, self.nTerms):
            comboterm.remove_text(0)
          self.nTerms = 0

          date = widget.get_active_text()
          if date is None:
              self.abort()
              return

          # Fill in the term's menu
          mData = self.fdb.getMetaData()
          arr = mData[self.output][date]
          numOldTerm = 0
          if len(arr) == 1:
            comboterm.append_text("/")
            self.nTerms += 1
          else:
            for t in arr:
              comboterm.append_text(t)
              if t == oldTerm:
               numOldTerm = self.nTerms
              self.nTerms += 1
          self.readyForUpdate = True
          comboterm.set_active(numOldTerm)

       except ValueError, Argument:
          self.abort()
          print "Error fired here !!", Argument

    #===========================================================================
    # onSelectTerm:
    # handles the event when a term is chosen, updates the choice of
    # the background and orders repainting plots.
    #===========================================================================
    def onSelectTerm(self, widget):
      if not self.readyForUpdate:
         return
      self.expe = None
      try:
        self.readyForUpdate = False
        comboexpe = self.wTree.get_widget("comboexpe")
        combodate = self.wTree.get_widget("combodate")
        if widget.get_active_text() is None:
            self.abort()
            return
        if comboexpe.get_active_text() is None:
            self.abort()
            return
        if combodate.get_active_text() is None:
            self.abort()
            return

        self.termTitle = ''
        if widget.get_active_text() == "/":
          expe = self.output + "." + combodate.get_active_text()
        else:
          self.termTitle = widget.get_active_text()
          expe = self.output + "." + combodate.get_active_text() + "." + widget.get_active_text()
        self.expeData = self.fdb.map[expe]

        bg = self.experiments.getBackground(self.idExpe)
        self.background = self.fdb.getExpeAt(bg['expe'], bg['refTime'], self.expeData["validTime"])
        self.readyForUpdate = True
        self.update()

      except ValueError, Argument:
        self.abort()
        print "Error fired here !!", Argument


    #===========================================================================
    # abort:
    # Abort current selection
    #===========================================================================          
    def abort(self):
      self.readyForUpdate = True
      self.expeData = None
      self.update()

    #===========================================================================
    # update:
    # Updates datas in the plot according to the choices done in the GUI.
    #===========================================================================
    def update(self):
      self.plt.clear()
      self.plot()
      self.canvas.draw()

    #===========================================================================
    # plot:
    # Plot data on the canvas
    #===========================================================================
    def plot(self):
      if self.expeData is  None:
        return

      obsData = None
      if self.odb is not None:
         obsData = self.odb.getXYAt(self.expeData["validTime"], "ObsVal")
      expeData = self.fdb.getData(self.expeData)

      bgData = None
      if self.background is not None:
         bgData = self.fdb.getData(self.background)

      title = "Experiment: " + self.idExpe + " at " + self.expeData["referenceTime"] + ' ' + self.termTitle
      self.plt.setPlot(title, expeData, bgData, obsData)
