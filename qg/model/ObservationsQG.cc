/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ObservationsQG.h"

#include <map>
#include <string>

#include "util/Logger.h"
#include "util/abor1_cpp.h"
#include "eckit/config/Configuration.h"


using oops::Log;

namespace qg {

// -----------------------------------------------------------------------------

// We have to use a pointer here because we don't know the order in
// which objects are created and this causes problems with xlc.

std::map < std::string, ObsFactory * > * ObsFactory::makers_ = 0;

// -----------------------------------------------------------------------------

ObsFactory::ObsFactory(const std::string & name) {
  if (!makers_) makers_=new std::map < std::string, ObsFactory * >();

  if (makers_->find(name) != makers_->end()) {
    Log::error() << name << " already registered in QG observation factory." << std::endl;
    ABORT("Element already registered in ObsFactory.");
  }
  (*makers_)[name] = this;
}

// -----------------------------------------------------------------------------

ObservationsQG * ObsFactory::create(ObsSpaceQG & odb, const eckit::Configuration & conf) {
  if (!makers_) makers_=new std::map < std::string, ObsFactory * >();

  std::string id = conf.getString("ObsType");
  Log::trace() << "QG ObsFactory ObsType =" << id << std::endl;
  std::map<std::string, ObsFactory*>::iterator j = makers_->find(id);
  if (j == makers_->end()) {
    Log::error() << id << " does not exist in QG observation factory." << std::endl;
    ABORT("Element does not exist in ObsFactory.");
  }
  return (*j).second->make(odb, conf);
}

// -----------------------------------------------------------------------------
}  // namespace qg
