/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#ifndef OOPS_BASE_INFLATIONBASE_H_
#define OOPS_BASE_INFLATIONBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/StateSet.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

class InflationParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(InflationParameters, Parameters);
 public:
    RequiredParameter<std::string> method{"method", this};
};

template<typename MODEL>
class InflationBase {
    typedef Geometry<MODEL>                   Geometry_;
    typedef IncrementSet<MODEL>               IncrementSet_;
    typedef StateSet<MODEL>                   StateSet_;
    typedef InflationParameters               Parameters_;

 public:
  InflationBase(const Geometry_ &,
                const StateSet_ &, const Variables &);
  virtual ~InflationBase() = default;

  virtual void doInflation(IncrementSet_ &) = 0;
  virtual void doInflation(StateSet_ &) = 0;

 protected:
  const Geometry_ & geometry() const {return geom_;}
  const StateSet_ & background() const {return bgEns_;}
  const Variables & vars() const {return anVars_;}

 private:
  const Geometry_ & geom_;
  const StateSet_ & bgEns_;
  const Variables & anVars_;
};


template<class MODEL>
class InflationFactory {
  typedef Geometry<MODEL>                   Geometry_;
  typedef StateSet<MODEL>                   StateSet_;

 public:
  static InflationBase<MODEL> * create(const eckit::LocalConfiguration &, const Geometry_ &,
                                const StateSet_ &, const Variables &);
  virtual ~InflationFactory() = default;

 protected:
  /// \brief Register a maker able to create observation operators of type \p name.
  explicit InflationFactory(const std::string &name);

 private:
  virtual InflationBase<MODEL> * make(const eckit::LocalConfiguration &, const Geometry_ &,
                                 const StateSet_ &, const Variables &) = 0;

  static std::map < std::string, InflationFactory * > & getMakers() {
    static std::map < std::string, InflationFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template <typename MODEL>
InflationFactory<MODEL>::InflationFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in inflation factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
InflationBase<MODEL>* InflationFactory<MODEL>::create(const eckit::LocalConfiguration & config,
                                                      const Geometry_ & geom,
                                                      const StateSet_ & bens,
                                                      const Variables & vars) {
  std::string id = config.getString("method");
  Log::trace() << "Inflation Type= " << id << std::endl;
  typename std::map<std::string, InflationFactory<MODEL>*>::iterator j = getMakers().find(id);
  if (j == getMakers().end()) {
    throw std::runtime_error(id + " does not exist in inflation factory.");
  }
  Log::trace() << "InflationFactory::create inflation type" << std::endl;
  return (*j).second->make(config, geom, bens, vars);
}

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class InflationMaker : public InflationFactory<MODEL> {
  typedef typename T::Parameters_ Parameters_;

  InflationBase<MODEL> * make(const eckit::LocalConfiguration & config,
                              const Geometry<MODEL> & geom,
                              const StateSet<MODEL> & bens,
                              const Variables & vars) override {
    Parameters_ params;
    params.deserialize(config);
    return new T(params, geom, bens, vars);
  }

 public:
  explicit InflationMaker(const std::string & name) : InflationFactory<MODEL>(name) {}
};

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

template <typename MODEL>
InflationBase<MODEL>::InflationBase(const Geometry_ & geom,
                                    const StateSet_ & bens,
                                    const Variables & vars)
  : geom_(geom), bgEns_(bens), anVars_(vars)
{}

}  // namespace oops
#endif  // OOPS_BASE_INFLATIONBASE_H_
