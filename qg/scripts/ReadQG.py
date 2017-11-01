# (C) Copyright 2009-2016 ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

import numpy as np

class QgRead:

  def readFields(self,fdbdir,attributes):
    filename = fdbdir + "/" +                                           \
           attributes["expver"]        +  "." + \
           attributes["type"]          +  "." + \
           attributes["referencedate"]

    if attributes["type"]=="fc": filename += "." + attributes["step"]

    print 'reading from file ' + filename
    infile = open(filename,"r")

    [nx, ny, nl, nf, ns] = [int(i) for i in infile.readline().split()]
    print "grid dimensions are: nx=",nx," ny=",ny

    validity_date = infile.readline()
    print "validity date is: ",validity_date

    psi = np.fromfile(file=infile,dtype="float64", \
                      count=(2*ny*nx),sep=" ").reshape((2,ny,nx))
    print "max(psi)=",np.amax(psi)," min(psi)=",np.amin(psi)

    pv  = np.fromfile(file=infile,dtype="float64", \
                      count=(2*ny*nx),sep=" ").reshape((2,ny,nx))
    print "max(pv)=",np.amax(pv)," min(pv)=",np.amin(pv)

    u   = np.fromfile(file=infile,dtype="float64",\
                      count=(2*ny*nx),sep=" ").reshape((2,ny,nx))
    print "max(u)=",np.amax(u)," min(u)=",np.amin(u)

    v   = np.fromfile(file=infile,dtype="float64", \
                      count=(2*ny*nx),sep=" ").reshape((2,ny,nx))
    print "max(v)=",np.amax(v)," min(v)=",np.amin(v)

    infile.close()

    return pv, psi
