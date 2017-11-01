/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_POSTPROCESSORTL_H_
#define OOPS_BASE_POSTPROCESSORTL_H_

#include <vector>

#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

#include "oops/base/PostBaseTL.h"
#include "oops/base/GeneralizedDepartures.h"

namespace oops {

/// Control model post processing
/*!
 *  This class controls model post processing in the most general sense,
 *  ie all diagnostics computations that do not affect the model integration.
 *  It just calls all the individual processors one by one.
 */

template<typename INCR>
class PostProcessorTL {
  typedef PostBaseTL<INCR> PostBaseTL_;

 public:
  PostProcessorTL() {}
  PostProcessorTL(const PostProcessorTL & pp): processors_(pp.processors_) {}
  ~PostProcessorTL() {}

  void enrollProcessor(PostBaseTL_ * pp) {
    if (pp != 0) {
      boost::shared_ptr<PostBaseTL_> sp(pp);
      processors_.push_back(sp);
    }
  }

  void enrollProcessor(boost::shared_ptr<PostBaseTL_> pp) {
    if (pp != 0) processors_.push_back(pp);
  }

  GeneralizedDepartures * releaseOutputFromTL(unsigned int ii) {
    GeneralizedDepartures * lambda = processors_[ii]->releaseOutputFromTL();
    return lambda;
  }

  void initializeTL(const INCR & dx, const util::DateTime & end,
                    const util::Duration & step) {
    BOOST_FOREACH(boost::shared_ptr<PostBaseTL_> jp, processors_) {
      jp->initializeTL(dx, end, step);
    }
  }

  void processTL(const INCR & dx) {
    BOOST_FOREACH(boost::shared_ptr<PostBaseTL_> jp, processors_) jp->processTL(dx);
  }

  void finalizeTL(const INCR & dx) {
    BOOST_FOREACH(boost::shared_ptr<PostBaseTL_> jp, processors_) jp->finalizeTL(dx);
  }

 private:
  std::vector< boost::shared_ptr<PostBaseTL_> > processors_;
  PostProcessorTL operator= (const PostProcessorTL &);
};

}  // namespace oops

#endif  // OOPS_BASE_POSTPROCESSORTL_H_
