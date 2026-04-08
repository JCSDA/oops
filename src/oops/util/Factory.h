/*
 * (C) Crown Copyright 2025-2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <string>

#include "oops/util/abor1_cpp.h"

namespace util {

// Forward declaration.
template<typename BaseClass, typename DerivedClass, typename... MakerArgs>
class FactoryMaker;

/// \brief Generic Factory design pattern.
/// \details Use with `using` to define a custom factory. For example:
///
///          ```
///          // Define factory.
///          using StateFactory = Factory<State,  // Base class.
///                                       // Constructor/Maker arguments.
///                                       const Geometry&,
///                                       const Variables&>;
///          ```
///
///          Then, in the derived class, `.cc` file, at the top, register a derived type
///          to the factory, with a given key string:
///
///          ```
///          static StateFactory::Maker<MyDerivedState> makerMyDerivedState_("my-derived");
///          ```
///
///          To then instantiate a given derived type when needed:
///
///          ```
///          const std::unique_ptr<State> myState(StateFactory::create("my-derived",
///                                                                    geom,    // Constructor args
///                                                                    vars));
///          ```
template<typename BaseClass, typename... MakerArgs>
class Factory {
 private:
  virtual BaseClass * make(MakerArgs... args) = 0;

  static std::map<std::string, Factory *> & getMakers() {
    static std::map<std::string, Factory *> makers_;
    return makers_;
  }

 public:
  template<typename Derived>
  using Maker = FactoryMaker<BaseClass, Derived, MakerArgs...>;

  static BaseClass * create(const std::string & name, MakerArgs... args) {
    auto j = getMakers().find(name);
    if (j == getMakers().end()) {
      ABORT(name + " does not exist in factory: "
            + std::string(typeid(Factory).name()) + ".");
    }
    return (*j).second->make(args...);
  }

  virtual ~Factory() = default;

 protected:
  /// \brief Register a maker
  explicit Factory(const std::string & name) {
    if (getMakers().find(name) != getMakers().end()) {
      ABORT(name + " already registered in factory: "
            + std::string(typeid(Factory).name()) + ".");
    }
    getMakers()[name] = this;
  }
};

// -----------------------------------------------------------------------------

template<typename BaseClass, typename DerivedClass, typename... MakerArgs>
class FactoryMaker : public Factory<BaseClass, MakerArgs...> {
  BaseClass * make(MakerArgs... args) override {
    return new DerivedClass(args...);
  }

 public:
  explicit FactoryMaker(const std::string & name) : Factory<BaseClass, MakerArgs...>(name) {}
};

}  // namespace util
