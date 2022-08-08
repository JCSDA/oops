/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/ObsOpBaseTLAD.h"

#include "eckit/config/Configuration.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace qg {

// -----------------------------------------------------------------------------

ObsOpTLADFactory::ObsOpTLADFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in qg::ObsOpTLADFactory." << std::endl;
    ABORT("Element already registered in ObsOpTLADFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

ObsOpBaseTLAD * ObsOpTLADFactory::create(const ObsSpaceQG & odb,
                                         const eckit::Configuration & conf) {
  oops::Log::trace() << "ObsOpBaseQG::create starting" << std::endl;
  const std::string id = conf.getString("obs type");
  typename std::map<std::string, ObsOpTLADFactory*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    oops::Log::error() << id << " does not exist in observation operator factory." << std::endl;
    ABORT("Element does not exist in ObsOpTLADFactory.");
  }
  ObsOpBaseTLAD * ptr = jloc->second->make(odb, conf);
  oops::Log::trace() << "ObsOpBaseQG::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace qg
