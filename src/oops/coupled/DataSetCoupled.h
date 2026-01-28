/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <vector>

#include "oops/base/DataSetBase.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment4D.h"
#include "oops/base/State4D.h"
#include "oops/coupled/IncrementCoupled.h"
#include "oops/coupled/StateCoupled.h"
#include "oops/coupled/TraitCoupled.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------
// Specialized IncrementSet constructor implementations for zero-copy component access

template<typename MODEL1, typename MODEL2>
Increment4D<MODEL1> share_increment1(const Increment4D<TraitCoupled<MODEL1, MODEL2>> & coupled) {
  std::vector<std::shared_ptr<Increment<MODEL1>>> incset;
  for (size_t jt = 0; jt < coupled.time_size(); ++jt) {
    incset.emplace_back(coupled[jt].increment().getIncrementPtr1());
  }
  Increment4D<MODEL1> inc(coupled.times(), coupled.commTime(),
                          coupled.members(), coupled.commEns(), incset);
  return inc;
}

template<typename MODEL1, typename MODEL2>
Increment4D<MODEL2> share_increment2(const Increment4D<TraitCoupled<MODEL1, MODEL2>> & coupled) {
  std::vector<std::shared_ptr<Increment<MODEL2>>> incset;
  for (size_t jt = 0; jt < coupled.time_size(); ++jt) {
    incset.emplace_back(coupled[jt].increment().getIncrementPtr2());
  }
  Increment4D<MODEL2> inc(coupled.times(), coupled.commTime(),
                          coupled.members(), coupled.commEns(), incset);
  return inc;
}

template<typename MODEL1, typename MODEL2>
State4D<MODEL1> share_state1(const State4D<TraitCoupled<MODEL1, MODEL2>> & coupled) {
  std::vector<std::shared_ptr<State<MODEL1>>> stateset;
  for (size_t jt = 0; jt < coupled.time_size(); ++jt) {
    stateset.emplace_back(coupled[jt].state().getStatePtr1());
  }
  State4D<MODEL1> state(coupled.times(), coupled.commTime(),
                        coupled.members(), coupled.commEns(), stateset);
  return state;
}

template<typename MODEL1, typename MODEL2>
State4D<MODEL2> share_state2(const State4D<TraitCoupled<MODEL1, MODEL2>> & coupled) {
  std::vector<std::shared_ptr<State<MODEL2>>> stateset;
  for (size_t jt = 0; jt < coupled.time_size(); ++jt) {
    stateset.emplace_back(coupled[jt].state().getStatePtr2());
  }
  State4D<MODEL2> state(coupled.times(), coupled.commTime(),
                        coupled.members(), coupled.commEns(), stateset);
  return state;
}

// -----------------------------------------------------------------------------

}  // namespace oops
