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
#include "oops/assimilation/DRPBlockLanczosMinimizer.h"
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
#include "oops/assimilation/PFF.h"
#include "oops/assimilation/PLanczosMinimizer.h"
#include "oops/assimilation/RPCGMinimizer.h"
#include "oops/assimilation/RPLanczosMinimizer.h"
#include "oops/assimilation/SaddlePointMinimizer.h"


namespace oops {

template <typename MODEL, typename OBS> void instantiateMinFactory() {
  static MinMaker<MODEL, OBS, DRGMRESRMinimizer<MODEL, OBS> >     makerDRGMRESR_("DRGMRESR");
  static MinMaker<MODEL, OBS, DRIPCGMinimizer<MODEL, OBS> >       makerDRIPCG_("DRIPCG");
  static MinMaker<MODEL, OBS, GMRESRMinimizer<MODEL, OBS> >       makerGMRESR_("GMRESR");
  static MinMaker<MODEL, OBS, IPCGMinimizer<MODEL, OBS> >         makerIPCG_("IPCG");
  static MinMaker<MODEL, OBS, SaddlePointMinimizer<MODEL, OBS> >  makerSADDLE_("SaddlePoint");
  static MinMaker<MODEL, OBS, RPCGMinimizer<MODEL, OBS> >         makerRPCG_("RPCG");
  static MinMaker<MODEL, OBS, DRPCGMinimizer<MODEL, OBS> >        makerDRPCG_("DRPCG");
  static MinMaker<MODEL, OBS, DRPFOMMinimizer<MODEL, OBS> >       makerDRPFOM_("DRPFOM");
  static MinMaker<MODEL, OBS, LBGMRESRMinimizer<MODEL, OBS> >     makerBDRPCG_("LBGMRESR");
  static MinMaker<MODEL, OBS, DRPLanczosMinimizer<MODEL, OBS> >   makerDRPLanczos_("DRPLanczos");
  static MinMaker<MODEL, OBS, PCGMinimizer<MODEL, OBS> >          makerPCG_("PCG");
  static MinMaker<MODEL, OBS, PLanczosMinimizer<MODEL, OBS> >     makerPLanczos_("PLanczos");
  static MinMaker<MODEL, OBS, RPLanczosMinimizer<MODEL, OBS> >    makerRPLanczos_("RPLanczos");
  static MinMaker<MODEL, OBS, MINRESMinimizer<MODEL, OBS> >       makerMINRES_("MINRES");
  static MinMaker<MODEL, OBS, FGMRESMinimizer<MODEL, OBS> >       makerFGMRES_("FGMRES");
  static MinMaker<MODEL, OBS, DRPBlockLanczosMinimizer<MODEL, OBS> >
            makerBlockBLanczos_("DRPBlockLanczos");
  static MinMaker<MODEL, OBS, PFF<MODEL, OBS> >                   makerPFF_("PFF");
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_INSTANTIATEMINFACTORY_H_
