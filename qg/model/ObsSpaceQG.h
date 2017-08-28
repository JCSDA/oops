/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_OBSSPACEQG_H_
#define QG_MODEL_OBSSPACEQG_H_

#include <map>
#include <ostream>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "util/DateTime.h"
#include "util/Logger.h"
#include "util/Printable.h"

#include "model/ObsHelpQG.h"
#include "model/QgFortran.h"

namespace qg {
  class ObsVecQG;

/// Wrapper around ObsHelpQG, mostly to hide the factory

class ObsSpaceQG : public util::Printable {
 public:
  ObsSpaceQG(const eckit::Configuration &, const util::DateTime &, const util::DateTime &);
  ObsSpaceQG(const ObsSpaceQG &);
  ~ObsSpaceQG();

  void getdb(const std::string & col, int & keyData) const {
    helper_->getdb(obsname_, col, keyData);
  }
  void putdb(const std::string & col, const int & keyData) const {
    helper_->putdb(obsname_, col, keyData);
  }

  int locations(const util::DateTime & t1, const util::DateTime & t2) const {
    int key_locs;
    key_locs = helper_->locations(obsname_, t1, t2);
    return key_locs;
  }

  void generateDistribution(const eckit::Configuration & conf) {
    helper_->generateDistribution(conf, obsname_, winbgn_, winend_, nobs_);
  }

  void printJo(const ObsVecQG &, const ObsVecQG &);

  int nobs() const {return nobs_;}
  int nvin() const {return nvin_;}
  int nout() const {return nout_;}
  const std::string & obsname() const {return obsname_;}
  const util::DateTime & windowStart() const {return winbgn_;}
  const util::DateTime & windowEnd() const {return winend_;}
  const eckit::Configuration & config() const {return conf_;}

  int & toFortran() {return helper_->toFortran();}
  const int & toFortran() const {return helper_->toFortran();}

 private:
  void print(std::ostream &) const;
  ObsSpaceQG & operator= (const ObsSpaceQG &);
  std::string ref_;
  mutable ObsHelpQG * helper_;
  std::string obsname_;
  unsigned int nobs_;
  unsigned int nvin_;
  unsigned int nout_;
  const eckit::LocalConfiguration conf_;
  const util::DateTime winbgn_;
  const util::DateTime winend_;

  static std::map < std::string, int > theObsFileCount_;
};

}  // namespace qg

#endif  // QG_MODEL_OBSSPACEQG_H_
