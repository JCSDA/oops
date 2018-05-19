/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_OBSHELPQG_H_
#define QG_MODEL_OBSHELPQG_H_

#include <string>

#include <boost/noncopyable.hpp>

#include "model/QgFortran.h"
#include "oops/util/ObjectCounter.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace qg {

/// Observation Data Handler for QG Model

class ObsHelpQG : private boost::noncopyable,
                  private util::ObjectCounter<ObsHelpQG> {
 public:
  static const std::string classname() {return "qg::ObsHelpQG";}

  explicit ObsHelpQG(const eckit::Configuration &);
  ~ObsHelpQG();

  void getdb(const std::string &, const std::string &, int & keyOvec) const;
  void putdb(const std::string &, const std::string &, const int & keyOvec);

  F90locs locations(const std::string &, const util::DateTime &, const util::DateTime &) const;
  void generateDistribution(const eckit::Configuration &, const std::string &,
                            const util::DateTime &, const util::DateTime &, unsigned int &);
  int nobs(const std::string &) const;

  int & toFortran() {return keyHelp_;}
  const int & toFortran() const {return keyHelp_;}

 private:
  F90odb keyHelp_;
};

}  // namespace qg

#endif  // QG_MODEL_OBSHELPQG_H_
