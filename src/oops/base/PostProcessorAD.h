/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_POSTPROCESSORAD_H_
#define OOPS_BASE_POSTPROCESSORAD_H_

#include <vector>

#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

#include "oops/base/PostBaseAD.h"
#include "oops/base/GeneralizedDepartures.h"

namespace oops {

/// Control model post processing
/*!
 *  This class controls model post processing in the most general sense,
 *  ie all diagnostics computations that do not affect the model integration.
 *  It just calls all the individual processors one by one.
 */

template<typename INCR>
class PostProcessorAD {
  typedef PostBaseAD<INCR> PostBaseAD_;

 public:
  PostProcessorAD() {}
  PostProcessorAD(const PostProcessorAD & pp): processors_(pp.processors_) {}
  ~PostProcessorAD() {}

  void enrollProcessor(PostBaseAD_ * pp) {
    if (pp != 0) {
      boost::shared_ptr<PostBaseAD_> sp(pp);
      processors_.push_back(sp);
    }
  }

  void enrollProcessor(boost::shared_ptr<PostBaseAD_> pp) {
    if (pp != 0) processors_.push_back(pp);
  }

  void initializeAD(INCR & dx, const util::DateTime & bgn,
                    const util::Duration & step) {
    BOOST_FOREACH(boost::shared_ptr<PostBaseAD_> jp, processors_) {
      jp->initializeAD(dx, bgn, step);
    }
  }

  void processAD(INCR & dx) {
    BOOST_FOREACH(boost::shared_ptr<PostBaseAD_> jp, processors_) jp->processAD(dx);
  }

  void finalizeAD(INCR & dx) {
    BOOST_FOREACH(boost::shared_ptr<PostBaseAD_> jp, processors_) jp->finalizeAD(dx);
  }

 private:
  std::vector< boost::shared_ptr<PostBaseAD_> > processors_;
  PostProcessorAD operator= (const PostProcessorAD &);
};

}  // namespace oops

#endif  // OOPS_BASE_POSTPROCESSORAD_H_
