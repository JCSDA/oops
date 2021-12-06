/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_POSTPROCESSORTLAD_H_
#define OOPS_BASE_POSTPROCESSORTLAD_H_

#include <vector>

#include <memory>

#include "oops/base/Increment.h"
#include "oops/base/PostBaseTLAD.h"
#include "oops/base/State.h"

namespace oops {

/// Control model post processing
/*!
 *  This class controls model post processing in the most general sense,
 *  ie all diagnostics computations that do not affect the model integration.
 *  It just calls all the individual processors one by one.
 */

template<typename MODEL>
class PostProcessorTLAD {
  typedef Increment<MODEL>           Increment_;
  typedef PostBaseTLAD<MODEL>        PostBaseTLAD_;
  typedef State<MODEL>               State_;

 public:
  PostProcessorTLAD() {}
  PostProcessorTLAD(const PostProcessorTLAD & pp): processors_(pp.processors_) {}
  ~PostProcessorTLAD() {}

  void enrollProcessor(PostBaseTLAD_ * pp) {
    if (pp != 0) {
      std::shared_ptr<PostBaseTLAD_> sp(pp);
      processors_.push_back(sp);
    }
  }

  void enrollProcessor(std::shared_ptr<PostBaseTLAD_> pp) {
    if (pp != 0) processors_.push_back(pp);
  }

/// Set linearization state
  void initializeTraj(const State_ & xx, const util::DateTime & end,
                      const util::Duration & step) {
    for (auto & jp : processors_) {
      jp->initializeTraj(xx, end, step);
    }
  }

  void processTraj(const State_ & xx) {
    for (auto & jp : processors_) {
      jp->processTraj(xx);
    }
  }

  void finalizeTraj(const State_ & xx) {
    for (auto & jp : processors_) {
      jp->finalizeTraj(xx);
    }
  }

/// Tangent linear methods
  void initializeTL(const Increment_ & dx, const util::DateTime & end,
                    const util::Duration & step) {
    for (auto & jp : processors_) {
      jp->initializeTL(dx, end, step);
    }
  }

  void processTL(const Increment_ & dx) {
    for (auto & jp : processors_) {
      jp->processTL(dx);
    }
  }

  void finalizeTL(const Increment_ & dx) {
    for (auto & jp : processors_) {
      jp->finalizeTL(dx);
    }
  }

/// Adjoint methods
  void initializeAD(Increment_ & dx, const util::DateTime & bgn,
                    const util::Duration & step) {
    for (auto & jp : processors_) {
      jp->initializeAD(dx, bgn, step);
    }
  }

  void processAD(Increment_ & dx) {
    for (auto & jp : processors_) {
      jp->processAD(dx);
    }
  }

  void finalizeAD(Increment_ & dx) {
    for (auto & jp : processors_) {
      jp->finalizeAD(dx);
    }
  }

 private:
  std::vector< std::shared_ptr<PostBaseTLAD_> > processors_;
  PostProcessorTLAD operator= (const PostProcessorTLAD &);
};

}  // namespace oops

#endif  // OOPS_BASE_POSTPROCESSORTLAD_H_
