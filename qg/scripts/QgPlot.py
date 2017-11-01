# (C) Copyright 2009-2016 ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class QgPlot:

  def plot(self,pv,psi):
    if (pv.shape != psi.shape):
      raise "pv.shape != psi.shape"

    nz,ny,nx=pv.shape

    x   = np.arange(1,nx+1)
    y   = np.arange(1,ny+1)

    # colour map
    CM = cm.jet

    plt.figure()

    # ================= upper plot (k=0)
    ax = plt.subplot(2,1,1) # 2 rows, 1 column, first plot
    ax.set_autoscale_on(False)
    plt.axis([0, nx+1, 0, ny+1])
    plt.title("Upper Layer")

    CSPV  = plt.contourf(x,y, pv[0,:,:],20,cmap=CM)
    CSPSI = plt.contour (x,y,psi[0,:,:],10,colors='k')

    plt.clabel(CSPSI, fontsize=9, inline=1)
    CB = plt.colorbar(CSPV)

    # ================= lower plot (k=1)
    bx = plt.subplot(2,1,2)
    bx.set_autoscale_on(False)
    plt.axis([0, nx+1, 0, ny+1])
    plt.title("Lower Layer")

    CSPV  = plt.contourf(x,y, pv[1,:,:],20,cmap=CM)
    CSPSI = plt.contour (x,y,psi[1,:,:],10,colors='k')

    plt.clabel(CSPSI, fontsize=9, inline=1)
    CB = plt.colorbar(CSPV)

    plt.show()
