/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_LOCATIONSQG_H_
#define QG_MODEL_LOCATIONSQG_H_

#include <ostream>
#include <string>

#include "model/ObsSpaceQG.h"
#include "model/QgFortran.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace qg {

/// LocationsQG class to handle locations for QG model.

class LocationsQG : public util::Printable,
              private util::ObjectCounter<LocationsQG> {
 public:
  static const std::string classname() {return "qg::LocationsQG";}

  LocationsQG(const ObsSpaceQG & ot,
        const util::DateTime & t1, const util::DateTime & t2) {
    keyLoc_ = ot.locations(t1, t2);
  }

  ~LocationsQG() {qg_loc_delete_f90(keyLoc_);}

  int nobs() const {
    int nobs;
    qg_loc_nobs_f90(keyLoc_, nobs);
    return nobs;
  }

  int toFortran() const {return keyLoc_;}
 private:
  void print(std::ostream & os) const {
    os << "LocationsQG::print not implemented";
  }
  F90locs keyLoc_;
};

}  // namespace qg

#endif  // QG_MODEL_LOCATIONSQG_H_
