/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_ENSEMBLESCOLLECTION_H_
#define OOPS_BASE_ENSEMBLESCOLLECTION_H_

#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>

#include "oops/base/Ensemble.h"
#include "util/DateTime.h"

namespace oops {

// -----------------------------------------------------------------------------

/// This class handles access to the ensembles of perturbation with a
/// DateTime

template<typename MODEL> class EnsemblesCollection {
  typedef boost::unordered_map<util::DateTime,
                               boost::shared_ptr<Ensemble<MODEL> > > EnsemblesMap;

 public:
/// Constructor

/// Destructor
  virtual ~EnsemblesCollection() {
  }

/// Accessors
  boost::shared_ptr<Ensemble<MODEL> > operator[](const util::DateTime& d) const {
    return ensMap_.at(d);
  }

  void put(const util::DateTime& d, boost::shared_ptr<Ensemble<MODEL> > e) {
    ensMap_[d] = e;
  }

  static EnsemblesCollection<MODEL>& getInstance() {
    static boost::shared_ptr<EnsemblesCollection<MODEL> > instance_;
    if ( instance_.use_count() == 0 )
      instance_.reset(new EnsemblesCollection<MODEL>());
    return *instance_;
  }

 private:
  EnsemblesCollection() {
  }

  EnsemblesMap ensMap_;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_ENSEMBLESCOLLECTION_H_
