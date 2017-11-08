/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSSPACES_H_
#define OOPS_BASE_OBSSPACES_H_

#include <cstddef>
#include <ostream>
#include <map>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/interface/ObservationSpace.h"
#include "util/DateTime.h"
#include "util/Logger.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "util/Timer.h"

namespace oops {
  template <typename T>
  class Departures;

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsSpaces : public util::Printable,
                  private boost::noncopyable,
                  private util::ObjectCounter<ObsSpaces<MODEL> > {
  typedef Departures<MODEL>         Departures_;
  typedef ObservationSpace<MODEL>   ObsSpace_;

 public:
  static const std::string classname() {return "oops::ObsSpaces";}

  ObsSpaces(const eckit::Configuration &, const util::DateTime &, const util::DateTime &);
  ~ObsSpaces();

/// Access
  std::size_t size() const {return spaces_.size();}
  const ObsSpace_ & operator[](const std::size_t ii) const {return *spaces_.at(ii);} 
  std::size_t itype(const std::string & type) const {return types_.at(type);}

/// Assimilation window
  const util::DateTime & windowStart() const {return wbgn_;}
  const util::DateTime & windowEnd() const {return wend_;}

/// Other
  void printJo(const Departures_ &, const Departures_ &) const;  // To be changed

 private:
  void print(std::ostream &) const;
  std::vector<boost::shared_ptr<ObsSpace_> > spaces_;
  std::map<std::string, std::size_t> types_;
  const util::DateTime wbgn_;
  const util::DateTime wend_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsSpaces<MODEL>::ObsSpaces(const eckit::Configuration & conf,
                            const util::DateTime & bgn, const util::DateTime & end)
 : spaces_(0), types_(), wbgn_(bgn), wend_(end)
{
  std::vector<eckit::LocalConfiguration> obsconf;
  conf.get("ObsTypes", obsconf);
  for (std::size_t jj = 0; jj < obsconf.size(); ++jj) {
    Log::debug() << "ObsSpaces::ObsSpaces : conf " << obsconf[jj] << std::endl;
    const std::string otype = obsconf[jj].getString("ObsType");
    ASSERT(types_.count(otype) == 0);
    types_[otype] = jj;
    boost::shared_ptr<ObsSpace_> tmp(new ObsSpace_(obsconf[jj], bgn, end));
    spaces_.push_back(tmp);
//  Generate locations etc... if required
    if (obsconf[jj].has("Generate")) {
      const eckit::LocalConfiguration gconf(obsconf[jj], "Generate");
      spaces_[jj]->generateDistribution(gconf);
    }
  }
  ASSERT(spaces_.size() >0);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsSpaces<MODEL>::~ObsSpaces() {}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsSpaces<MODEL>::print(std::ostream & os) const {
  for (std::size_t jj = 0; jj < spaces_.size(); ++jj) {
    os << *spaces_[jj];
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsSpaces<MODEL>::printJo(const Departures_ & dy, const Departures_ & grad) const {
  for (std::size_t jj = 0; jj < spaces_.size(); ++jj) {
    spaces_[jj]->printJo(dy[jj], grad[jj]);
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSSPACES_H_
