/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_MODELDATA_H_
#define OOPS_INTERFACE_MODELDATA_H_

#include <memory>
#include <string>

#include "oops/base/Geometry.h"
#include "oops/base/Variables.h"

#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Model Data
///
template <typename MODEL>
class ModelData : public  util::Printable,
                  private util::ObjectCounter<ModelData<MODEL> >  {
  typedef typename MODEL::ModelData ModelData_;
  typedef Geometry<MODEL>           Geometry_;

 public:
  static const std::string classname() {return "oops::ModelData";}
  static const Variables defaultVariables() {return ModelData_::defaultVariables();}

  explicit ModelData(const Geometry_ &);
  virtual ~ModelData();
  ModelData(const ModelData &) = delete;
  ModelData(ModelData &&) = default;
  const ModelData & operator=(const ModelData &) = delete;
  ModelData & operator=(ModelData &&) = default;

  const eckit::LocalConfiguration modelData() const;

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<ModelData_> modeldata_;  // pointer to the ModelData implementation
};

// -----------------------------------------------------------------------------

template<typename MODEL>
ModelData<MODEL>::ModelData(const Geometry_ & geometry)
{
  Log::trace() << "ModelData<MODEL>::ModelData starting" << std::endl;
  util::Timer timer(classname(), "ModelData");
  modeldata_.reset(new ModelData_(geometry.geometry()));
  Log::trace() << "ModelData<MODEL>::ModelData done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ModelData<MODEL>::~ModelData() {
  Log::trace() << "ModelData<MODEL>::~ModelData starting" << std::endl;
  util::Timer timer(classname(), "~ModelData");
  modeldata_.reset();
  Log::trace() << "ModelData<MODEL>::~ModelData done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
const eckit::LocalConfiguration ModelData<MODEL>::modelData() const {
  Log::trace() << "ModelData<MODEL>::modelData starting" << std::endl;
  util::Timer timer(classname(), "modelData");
  Log::trace() << "ModelData<MODEL>::modelData done" << std::endl;
  return modeldata_->modelData();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ModelData<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ModelData<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *modeldata_;
  Log::trace() << "ModelData<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_MODELDATA_H_
