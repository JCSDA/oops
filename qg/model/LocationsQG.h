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

#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "model/QgFortran.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace qg {

// -----------------------------------------------------------------------------
/// LocationsQG class to handle locations for QG model.
class LocationsQG : public util::Printable,
                    private util::ObjectCounter<LocationsQG> {
 public:
  static const std::string classname() {return "qg::LocationsQG";}

// Constructors and basic operators
  explicit LocationsQG(const F90locs key) : keyLocs_(key) {}
  LocationsQG(const eckit::Configuration &, const eckit::mpi::Comm &);
  ~LocationsQG() {qg_locs_delete_f90(keyLocs_);}

// Utilities
  int size() const;
  int toFortran() const {return keyLocs_;}

 private:
  void print(std::ostream &) const;
  F90locs keyLocs_;
};

}  // namespace qg

#endif  // QG_MODEL_LOCATIONSQG_H_
