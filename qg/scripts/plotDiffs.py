#!/usr/local/bin/python
# (C) Copyright 2009-2016 ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


import QgPlot
import QgRead

#fdbdir="/home/rd/dai/Oops/Data/FDB"
fdbdir="/home/da/dasg/Data/FDB"
attributes = {}
attributes["expver"]        = "example"
attributes["validitydate"]  = "2010-01-01T00:00:00Z"
attributes["step"]  = "P1D"
qgReader = QgRead.QgRead()

# Read first fields
              
#background
attributes["type"]          = "fc"
attributes["referencedate"] = "2009-12-31T00:00:00Z"

#truth
#attributes["referencedate"] = "2009-12-15T00:00:00Z"
#attributes["type"]          = "truth"


print "reading first fields..."
pv1, psi1 = qgReader.readFields(fdbdir,attributes)

# Read second fields

#analysis
attributes["type"]          = "an"
attributes["referencedate"] = "2010-01-02T00:00:00Z"

#bg
#attributes["type"]          = "fc"
#attributes["referencedate"] = "2009-12-31T00:00:00Z"

print "reading second fields..."
pv2, psi2 = qgReader.readFields(fdbdir,attributes)

# Plot difference

qgplt = QgPlot.QgPlot()

qgplt.plot(pv2-pv1,psi2-psi1)
