/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_INFLATION_H_
#define OOPS_BASE_INFLATION_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/InflationBase.h"
#include "oops/base/instantiateInflationFactory.h"
#include "oops/base/State.h"
#include "oops/base/StateSet.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/ConfigHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace oops {

template<typename MODEL, class T>
class Inflation {
    typedef Geometry<MODEL>                   Geometry_;
    typedef IncrementSet<MODEL>               IncrementSet_;
    typedef StateSet<MODEL>                   StateSet_;

 public:
  Inflation(const eckit::LocalConfiguration &, const Geometry_ &,
            const StateSet_ &);
  Inflation(const eckit::LocalConfiguration &, const Geometry_ &,
            const StateSet_ &, const Variables &);

  void calculate(const std::vector<eckit::LocalConfiguration> &);

  void save(const eckit::Configuration &);

 private:
  const Geometry_ & geom_;
  const StateSet_ & bgens_;
  std::unique_ptr<T> anEns_;
};

// -----------------------------------------------------------------------------

template<typename MODEL, class T>
Inflation<MODEL, T>::Inflation(const eckit::LocalConfiguration & config, const Geometry_ & geom,
                               const StateSet_ & bgens)
    : geom_(geom), bgens_(bgens), anEns_()
{
  ASSERT((std::is_same_v<T, StateSet_>));
  anEns_.reset(new T(geom, config));
}

// -----------------------------------------------------------------------------

template<typename MODEL, class T>
Inflation<MODEL, T>::Inflation(const eckit::LocalConfiguration & config, const Geometry_ & geom,
                               const StateSet_ & bgens, const Variables & vars)
    : geom_(geom), bgens_(bgens), anEns_()
{
  ASSERT((std::is_same_v<T, IncrementSet_>));
  anEns_.reset(new T(geom, vars, bgens_.validTimes(), config, bgens_.commTime()));
}

// -----------------------------------------------------------------------------
template<typename MODEL, class T>
void Inflation<MODEL, T>::calculate(const std::vector<eckit::LocalConfiguration> & subconfigs) {
  for (size_t jj = 0; jj < subconfigs.size(); ++jj) {
    std::unique_ptr<InflationBase<MODEL>>
    infMethod(InflationFactory<MODEL>::create(subconfigs[jj], geom_, bgens_, anEns_->variables()));

    infMethod->doInflation(*anEns_);
    }
}

// -----------------------------------------------------------------------------

template<typename MODEL, class T>
void Inflation<MODEL, T>::save(const eckit::Configuration & config) {
  eckit::LocalConfiguration outConfig(config);
  for (size_t ii=0; ii < anEns_->size(); ++ii) {
    util::setMember(outConfig, ii+1);
    (*anEns_)[ii].write(outConfig);
  }
}

// -----------------------------------------------------------------------------


}  // namespace oops
#endif  // OOPS_BASE_INFLATION_H_
