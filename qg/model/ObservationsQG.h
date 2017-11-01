/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_OBSERVATIONSQG_H_
#define QG_MODEL_OBSERVATIONSQG_H_

#include <map>
#include <ostream>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "util/Printable.h"

namespace util {
  class DateTime;
}

namespace qg {
  class GomQG;
  class LinearObsOp;
  class ObsSpaceQG;
  class ObsVecQG;
  class ObsBias;
  class VariablesQG;

/// Observations for QG model.
/*!
 *  Base class for QG observations.
 */

// -----------------------------------------------------------------------------

class ObservationsQG;

/// ObsFactory
class ObsFactory {
 public:
  static ObservationsQG * create(ObsSpaceQG &, const eckit::Configuration &);
  virtual ~ObsFactory() { makers_->clear(); }
 protected:
  explicit ObsFactory(const std::string &);
 private:
  virtual ObservationsQG * make(ObsSpaceQG &, const eckit::Configuration &) =0;
  static std::map < std::string, ObsFactory * > * makers_;
};

template<class T>
class ObsMaker : public ObsFactory {
  virtual ObservationsQG * make(ObsSpaceQG & odb, const eckit::Configuration & c) {return new T(odb, c);}
 public:
  explicit ObsMaker(const std::string & name) : ObsFactory(name) {}
};

// -----------------------------------------------------------------------------

class ObservationsQG : public util::Printable,
                      private boost::noncopyable {
 public:
  static ObservationsQG * create(ObsSpaceQG & odb, const eckit::Configuration & conf)
    {return ObsFactory::create(odb, conf);}

  virtual ~ObservationsQG() {}

// Obs Operators
  virtual void obsEquiv(const GomQG &, ObsVecQG &, const ObsBias &) const =0;

// Get TLAD obs operator (should be more like static create(...)?
  virtual LinearObsOp * getTLAD() const =0;

// Other
  virtual void generateObsError(const eckit::Configuration &) =0;
  virtual boost::shared_ptr<const VariablesQG> variables() const =0;

  virtual int & toFortran() =0;
  virtual const int & toFortran() const =0;

 private:
  virtual void print(std::ostream &) const =0;
};

// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_OBSERVATIONSQG_H_
