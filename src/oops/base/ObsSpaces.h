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
#include <map>
#include <ostream>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/geometry/Point2.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

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
  ObsSpaces(const ObsSpaces &, const eckit::geometry::Point2 &, const double &, const int &);
  ~ObsSpaces();

/// Access
  std::size_t size() const {return spaces_.size();}
  const ObsSpace_ & operator[](const std::size_t ii) const {return *spaces_.at(ii);}

/// Assimilation window
  const util::DateTime & windowStart() const {return wbgn_;}
  const util::DateTime & windowEnd() const {return wend_;}

/// Other
  void printJo(const Departures_ &, const Departures_ &) const;  // To be changed

 private:
  void print(std::ostream &) const;
  std::vector<boost::shared_ptr<ObsSpace_> > spaces_;
  const util::DateTime wbgn_;
  const util::DateTime wend_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsSpaces<MODEL>::ObsSpaces(const eckit::Configuration & conf,
                            const util::DateTime & bgn, const util::DateTime & end)
  : spaces_(0), wbgn_(bgn), wend_(end)
{
  int member = conf.getInt("member", 0);
  std::vector<eckit::LocalConfiguration> typeconfs;
  conf.get("ObsTypes", typeconfs);
  for (std::size_t jj = 0; jj < typeconfs.size(); ++jj) {
    eckit::LocalConfiguration obsconf(typeconfs[jj], "ObsSpace");
    if (member) obsconf.set("member", member);
    Log::debug() << "ObsSpaces::ObsSpaces : conf " << obsconf << std::endl;
    boost::shared_ptr<ObsSpace_> tmp(new ObsSpace_(obsconf, bgn, end));
    spaces_.push_back(tmp);
//  Generate locations etc... if required
    if (typeconfs[jj].has("Generate")) {
      const eckit::LocalConfiguration gconf(typeconfs[jj], "Generate");
      spaces_[jj]->generateDistribution(gconf);
    }
  }
  ASSERT(spaces_.size() >0);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsSpaces<MODEL>::ObsSpaces(const ObsSpaces<MODEL> & obss, const eckit::geometry::Point2 & center,
                            const double & dist, const int & maxn)
  : spaces_(0), wbgn_(obss.wbgn_), wend_(obss.wend_)
{
  for (std::size_t jj = 0; jj < obss.size(); ++jj) {
    boost::shared_ptr<ObsSpace_> tmp(new ObsSpace_(obss[jj], center,
        dist, maxn));
    spaces_.push_back(tmp);
  }
  ASSERT(spaces_.size() == obss.size());
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
