/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSFILTERS_H_
#define OOPS_BASE_OBSFILTERS_H_

#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/FilterBase.h"
#include "oops/base/ObsFilter.h"
#include "oops/base/ObsSpaces.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsFilters : public util::Printable {
  typedef FilterBase<MODEL>         FilterBase_;
  typedef ObsFilter<MODEL>          ObsFilter_;
  typedef ObsSpaces<MODEL>          ObsSpace_;

 public:
  ObsFilters(const ObsSpace_ &, const eckit::Configuration &);
  ~ObsFilters();

  const ObsFilter_ & operator[](const std::size_t ii) const {return *filters_.at(ii);}

 private:
  void print(std::ostream &) const;
  std::vector<boost::shared_ptr<ObsFilter_> > filters_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsFilters<MODEL>::ObsFilters(const ObsSpace_ & os, const eckit::Configuration & conf)
  : filters_(os.size())
{
  for (std::size_t jj = 0; jj < os.size(); ++jj) {
    boost::shared_ptr<ObsFilter_> tmp(new ObsFilter_);
    filters_[jj] = tmp;
  }

  std::vector<eckit::LocalConfiguration> confs;
  conf.get("ObsFilters", confs);

  for (std::size_t jj = 0; jj < confs.size(); ++jj) {
    const std::string type = confs[jj].getString("ObsType");
    const std::size_t ii = os.itype(type);
    boost::shared_ptr<FilterBase_> tmp(FilterFactory<MODEL>::create(confs[jj]));
    filters_[ii]->enrollFilter(tmp);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsFilters<MODEL>::~ObsFilters() {}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsFilters<MODEL>::print(std::ostream & os) const {
  os << "ObsFilters for " << filters_.size() << " types:" << std::endl;
  for (std::size_t jj = 0; jj < filters_.size(); ++jj) {
    os << *filters_[jj] << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSFILTERS_H_
