/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_POSTPROCESSOR_H_
#define OOPS_BASE_POSTPROCESSOR_H_

#include <vector>

#include <memory>

#include "oops/base/PostBase.h"

namespace oops {

/// Control model post processing
/*!
 *  This class controls model post processing in the most general sense,
 *  ie all diagnostics computations that do not affect the model integration.
 *  It just calls all the individual processors one by one.
 */

template<typename FLDS>
class PostProcessor {
  typedef PostBase<FLDS> PostBase_;

 public:
  PostProcessor() {}
  PostProcessor(const PostProcessor & pp): processors_(pp.processors_) {}
  ~PostProcessor() {}

  void enrollProcessor(PostBase_ * pp) {
    if (pp != 0) {
      std::shared_ptr<PostBase_> sp(pp);
      processors_.push_back(sp);
    }
  }

  void enrollProcessor(std::shared_ptr<PostBase_> pp) {
    if (pp != 0) processors_.push_back(pp);
  }

  void initialize(const FLDS & xx, const util::DateTime & end,
                  const util::Duration & step) {
    for (auto & jp : processors_) {
      jp->initialize(xx, end, step);
    }
  }

  void process(const FLDS & xx) {
    for (auto & jp : processors_) {
      jp->process(xx);
    }
  }

  void finalize(const FLDS & xx) {
    for (auto & jp : processors_) {
      jp->finalize(xx);
    }
  }

 private:
  std::vector< std::shared_ptr<PostBase_> > processors_;
  PostProcessor operator= (const PostProcessor &);
};

}  // namespace oops

#endif  // OOPS_BASE_POSTPROCESSOR_H_
