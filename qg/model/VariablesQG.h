/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QG_MODEL_VARIABLESQG_H_
#define QG_MODEL_VARIABLESQG_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace qg {

// -----------------------------------------------------------------------------

class VariablesQG : public util::Printable,
                    private util::ObjectCounter<VariablesQG> {
 public:
  static const std::string classname() {return "qg::VariablesQG";}

  explicit VariablesQG(const oops::Variables &);
  explicit VariablesQG(const eckit::Configuration &);

  ~VariablesQG();

  VariablesQG(const VariablesQG &);

  const int * toFortran() const {return &fvars_[0];}
  const oops::Variables& toOopsVariables() const {return oopsvars_;}
 private:
  void print(std::ostream &) const;
  void setF90(const std::vector<std::string>);
  std::vector<int> fvars_;
  const oops::Variables oopsvars_;
};

// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_VARIABLESQG_H_
