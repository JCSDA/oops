/*
 * (C) Copyright 2020- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_INTERPOLATORATLAS_H_
#define OOPS_GENERIC_INTERPOLATORATLAS_H_

#include <memory>
#include <ostream>
#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/interpolation.h"
#include "eckit/config/Configuration.h"

#include "oops/base/InterpolatorBase.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/ObjectCounter.h"

namespace oops {

// -----------------------------------------------------------------------------

/*! \brief Interface for Atlas interpolation
 *
 */

class InterpolatorAtlas : public InterpolatorBase,
                          private util::ObjectCounter<InterpolatorAtlas> {
 public:
  static const std::string classname() {return "oops::InterpolatorAtlas";}

  InterpolatorAtlas(const eckit::Configuration &, const atlas::FunctionSpace &,
                    const atlas::FunctionSpace &,
                    const atlas::field::FieldSetImpl * = nullptr);

  ~InterpolatorAtlas() { }
  void apply(atlas::Field const &, atlas::Field &);
  void apply(atlas::FieldSet const &, atlas::FieldSet &) override;

  void apply_ad(atlas::FieldSet const &, atlas::FieldSet &) override {
    std::string ErrorMsg = "Interpolator Adjoint not yet implemented for Atlas";
    ABORT(ErrorMsg);
  }

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<atlas::Interpolation> interpolator_;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_INTERPOLATORATLAS_H_
