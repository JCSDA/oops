/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/WriteParametersBase.h"

namespace oops {

void WriteParametersBase::setMember(int mem) {
  eckit::LocalConfiguration conf = toConfiguration();
  conf.set("member", mem);
  validateAndDeserialize(conf);
}

void WriteParametersBase::setDate(const util::DateTime &d) {
  eckit::LocalConfiguration conf = toConfiguration();
  conf.set("date", d.toString());
  validateAndDeserialize(conf);
}

void WriteParametersBase::setIteration(int it) {
  eckit::LocalConfiguration conf = toConfiguration();
  conf.set("iteration", it);
  validateAndDeserialize(conf);
}

}  // namespace oops
