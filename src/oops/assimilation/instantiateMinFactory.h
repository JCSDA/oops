/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_INSTANTIATEMINFACTORY_H_
#define OOPS_ASSIMILATION_INSTANTIATEMINFACTORY_H_

#include "oops/assimilation/DRGMRESRMinimizer.h"
#include "oops/assimilation/DRIPCGMinimizer.h"
#include "oops/assimilation/DRPCGMinimizer.h"
#include "oops/assimilation/DRPFOMMinimizer.h"
#include "oops/assimilation/DRPLanczosMinimizer.h"
#include "oops/assimilation/FGMRESMinimizer.h"
#include "oops/assimilation/GMRESRMinimizer.h"
#include "oops/assimilation/IPCGMinimizer.h"
#include "oops/assimilation/LBGMRESRMinimizer.h"
#include "oops/assimilation/Minimizer.h"
#include "oops/assimilation/MINRESMinimizer.h"
#include "oops/assimilation/PCGMinimizer.h"
#include "oops/assimilation/PLanczosMinimizer.h"
#include "oops/assimilation/RPCGMinimizer.h"
#include "oops/assimilation/RPLanczosMinimizer.h"
#include "oops/assimilation/SaddlePointMinimizer.h"


namespace oops {

template <typename MODEL> void instantiateMinFactory() {
  static MinMaker<MODEL, DRGMRESRMinimizer<MODEL> >     makerDRGMRESR_("DRGMRESR");
  static MinMaker<MODEL, DRIPCGMinimizer<MODEL> >       makerDRIPCG_("DRIPCG");
  static MinMaker<MODEL, GMRESRMinimizer<MODEL> >       makerGMRESR_("GMRESR");
  static MinMaker<MODEL, IPCGMinimizer<MODEL> >         makerIPCG_("IPCG");
  static MinMaker<MODEL, SaddlePointMinimizer<MODEL> >  makerSADDLE_("SaddlePoint");
  static MinMaker<MODEL, RPCGMinimizer<MODEL> >         makerRPCG_("RPCG");
  static MinMaker<MODEL, DRPCGMinimizer<MODEL> >        makerDRPCG_("DRPCG");
  static MinMaker<MODEL, DRPFOMMinimizer<MODEL> >       makerDRPFOM_("DRPFOM");
  static MinMaker<MODEL, LBGMRESRMinimizer<MODEL> >     makerBDRPCG_("LBGMRESR");
  static MinMaker<MODEL, DRPLanczosMinimizer<MODEL> >   makerDRPLanczos_("DRPLanczos");
  static MinMaker<MODEL, PCGMinimizer<MODEL> >          makerPCG_("PCG");
  static MinMaker<MODEL, PLanczosMinimizer<MODEL> >     makerPLanczos_("PLanczos");
  static MinMaker<MODEL, RPLanczosMinimizer<MODEL> >    makerRPLanczos_("RPLanczos");
  static MinMaker<MODEL, MINRESMinimizer<MODEL> >       makerMINRES_("MINRES");
  static MinMaker<MODEL, FGMRESMinimizer<MODEL> >       makerFGMRES_("FGMRES");
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_INSTANTIATEMINFACTORY_H_
