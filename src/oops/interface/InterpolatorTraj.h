/*
 * (C) Copyright 2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_INTERFACE_INTERPOLATORTRAJ_H_
#define OOPS_INTERFACE_INTERPOLATORTRAJ_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class InterpolatorTraj : public util::Printable,
                         private util::ObjectCounter<InterpolatorTraj<MODEL> > {
  typedef typename MODEL::InterpolatorTraj      InterpolatorTraj_;

 public:
  static const std::string classname() {return "oops::InterpolatorTraj";}

  InterpolatorTraj();
  ~InterpolatorTraj();

/// Interfacing
  InterpolatorTraj_ & interpolatortraj() {return *traj_;}
  const InterpolatorTraj_ & interpolatortraj() const {return *traj_;}

 private:
  InterpolatorTraj(const InterpolatorTraj &);
  InterpolatorTraj & operator=(const InterpolatorTraj &);
  void print(std::ostream &) const;
  boost::scoped_ptr<InterpolatorTraj_> traj_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
InterpolatorTraj<MODEL>::InterpolatorTraj(): traj_() {
  Log::trace() << "InterpolatorTraj<MODEL>::InterpolatorTraj starting" << std::endl;
  util::Timer timer(classname(), "InterpolatorTraj");
  traj_.reset(new InterpolatorTraj_());
  Log::trace() << "InterpolatorTraj<MODEL>::InterpolatorTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
InterpolatorTraj<MODEL>::~InterpolatorTraj() {
  Log::trace() << "InterpolatorTraj<MODEL>::~InterpolatorTraj starting" << std::endl;
  util::Timer timer(classname(), "~InterpolatorTraj");
  traj_.reset();
  Log::trace() << "InterpolatorTraj<MODEL>::~InterpolatorTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void InterpolatorTraj<MODEL>::print(std::ostream & os) const {
  Log::trace() << "InterpolatorTraj<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *traj_;
  Log::trace() << "InterpolatorTraj<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_INTERPOLATORTRAJ_H_
