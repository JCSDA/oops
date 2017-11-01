/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "oops/base/WeightingFct.h"

#include <map>
#include <string>

#include "util/Logger.h"
#include "util/abor1_cpp.h"
#include "eckit/config/Configuration.h"

namespace oops {

// -----------------------------------------------------------------------------
//  Factory
// -----------------------------------------------------------------------------

// We have to use a pointer here because we don't know the order in
// which objects are created and this causes problems with xlc.

std::map < std::string, WeightFactory * > * WeightFactory::makers_ = 0;

// -----------------------------------------------------------------------------

WeightFactory::WeightFactory(const std::string & name) {
  if (!makers_) makers_=new std::map < std::string, WeightFactory * >();

  if (makers_->find(name) != makers_->end()) {
    Log::error() << name << " already registered in weight function factory." << std::endl;
    ABORT("Element already registered in WeightFactory.");
  }
  (*makers_)[name] = this;
}

// -----------------------------------------------------------------------------

WeightingFct * WeightFactory::create(const eckit::Configuration & config) {
  if (!makers_) makers_=new std::map < std::string, WeightFactory * >();

  std::string id = config.getString("type");
  std::map<std::string, WeightFactory *>::iterator j = makers_->find(id);
  if (j == makers_->end()) {
    Log::error() << id << " does not exist in weight function factory." << std::endl;
    ABORT("Element does not exist in WeightFactory.");
  }
  return (*j).second->make(config);
}

// -----------------------------------------------------------------------------

}  // namespace oops
