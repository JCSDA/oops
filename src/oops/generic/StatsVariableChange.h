/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_STATVARIABLECHANGE_H_
#define OOPS_GENERIC_STATVARIABLECHANGE_H_

#include <string>
#include <boost/noncopyable.hpp>
#include <vector>
#include <string>
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

/// Derived class of generic variable transform for statistical

template <typename MODEL>
class StatsVariableChange : public VariableChangeBase<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "oops::StatsVariableChange";}

  explicit StatsVariableChange(const eckit::Configuration &);
  virtual ~StatsVariableChange();

  void linearize(const State_ &,  const Geometry_ &) override;
  void multiply(const Increment_ &, Increment_ &) const override;
  void multiplyInverse(const Increment_ &, Increment_ &) const override;
  void multiplyAD(const Increment_ &, Increment_ &) const override;
  void multiplyInverseAD(const Increment_ &, Increment_ &) const override;

private:
  const std::vector<std::string> pathStrList;

  void print(std::ostream &) const override;
  std::string ExtractFilename(std::string const pathString);
  void ExtractModelVarForCalc(std::string const fileString,
                              std::string& modelVarToCalcString,
                              std::string& varRegrByString);
  void VariableLists(const std::vector<std::string> pathStrList, 
                     std::vector<std::string>& modelVarToCalcList, 
                     std::vector<std::string>& varRegrByList);

  //StatsVarData populate(const varin_ , const varout_, const eckit::Configuration);
  
};

template<typename MODEL>
StatsVariableChange<MODEL>::StatsVariableChange(const eckit::Configuration & conf)
  : VariableChangeBase<MODEL>(conf)
{
  Log::trace() << "StatsVariableChange<MODEL>::StatsVariableChange starting" << std::endl;
  Log::trace() << "StatsVariableChange<MODEL>::StatsVariableChange done" << std::endl;
}

template<typename MODEL>
StatsVariableChange<MODEL>::~StatsVariableChange() {
  Log::trace() << "StatsVariableChange<MODEL>::~StatsVariableChange starting" << std::endl;
  Log::trace() << "StatsVariableChange<MODEL>::~StatsVariableChange done" << std::endl;
}

template<typename MODEL>
void StatsVariableChange<MODEL>::linearize(const State_ & x,  const Geometry_ & resol) {
  Log::trace() << "StatsVariableChange<MODEL>::linearize starting" << std::endl;
  Log::trace() << "StatsVariableChange<MODEL>::linearize done" << std::endl;
}

template<typename MODEL>
void StatsVariableChange<MODEL>::multiply(const Increment_ & in, Increment_ & out) const {
  Log::trace() << "StatsVariableChange<MODEL>::multiply starting" << std::endl;

  UnstructuredGrid ug;
  in.convert_to(ug);
  out.convert_from(ug);

  Log::trace() << "StatsVariableChange<MODEL>::multiply done" << std::endl;
}

template<typename MODEL>
void StatsVariableChange<MODEL>::multiplyInverse(const Increment_ & in, Increment_ & out) const {
  Log::trace() << "StatsVariableChange<MODEL>::multiplyInverse starting" << std::endl;

  UnstructuredGrid ug;
  in.convert_to(ug);
  out.convert_from(ug);

  Log::trace() << "StatsVariableChange<MODEL>::multiplyInverse done" << std::endl;
}

template<typename MODEL>
void StatsVariableChange<MODEL>::multiplyAD(const Increment_ & in, Increment_ & out) const {
  Log::trace() << "StatsVariableChange<MODEL>::multiplyAD starting" << std::endl;

  UnstructuredGrid ug;
  in.convert_to(ug);
  out.convert_from(ug);

  Log::trace() << "StatsVariableChange<MODEL>::multiplyAD done" << std::endl;
}

template<typename MODEL>
void StatsVariableChange<MODEL>::multiplyInverseAD(const Increment_ & in, Increment_ & out) const {
  Log::trace() << "StatsVariableChange<MODEL>::multiplyInverseAD starting" << std::endl;

  UnstructuredGrid ug;
  in.convert_to(ug);
  out.convert_from(ug);

  Log::trace() << "StatsVariableChange<MODEL>::multiplyInverseAD done" << std::endl;
}

template<typename MODEL>
void StatsVariableChange<MODEL>::print(std::ostream & os) const {
}

template<typename MODEL>
std::string StatsVariableChange<MODEL>::ExtractFilename(std::string const pathString)
{
    int tempIndex;
    for (int iCharIndex = pathString.length(); iCharIndex > 0; --iCharIndex)
    {
        if (pathString.substr(iCharIndex - 1 ,1) == "/")
        {
             return pathString.substr(iCharIndex, pathString.length() - iCharIndex);
        }
        tempIndex = iCharIndex; 
    }
    std::cout << "ExtractFileName: tempIndex = " << tempIndex << std::endl;  
    throw std::runtime_error("FileName not extracted from path " + pathString);
} 

// extract {model_variable_to_be_calculated} string and {variable_regressed_by} string
// from fileString

template<typename MODEL>
void StatsVariableChange<MODEL>::ExtractModelVarForCalc(std::string const fileString, 
                                std::string& modelVarToCalcString, 
                                std::string& varRegrByString)
{
// FileName is of the form {modelVarToCalcString}__{varRegrByString}__reg.nc
// find first "__" and the last "__" to find the appropriate strings
    int firstUnderscoreIndex, lastUnderscoreIndex;
    for (int iCharIndex = 0;  iCharIndex < fileString.length() - 1; ++iCharIndex)
    {
        firstUnderscoreIndex = iCharIndex;
        if (fileString.substr(iCharIndex, 2) ==  "__")
        {
          break;
        }
    }
    for (int iCharIndex = fileString.length(); iCharIndex > 0; --iCharIndex)
    {
        lastUnderscoreIndex = iCharIndex - 1;
        if (fileString.substr(iCharIndex-1, 2) == "__")
        {
          break; 
        }
    }
    modelVarToCalcString = fileString.substr(0, firstUnderscoreIndex);
    varRegrByString = fileString.substr(firstUnderscoreIndex + 2, 
                                        lastUnderscoreIndex - firstUnderscoreIndex - 2);

}


// append the modelVarToCalcList and varRegrByList
template<typename MODEL>
void StatsVariableChange<MODEL>::VariableLists(const std::vector<std::string> pathStrList, 
                   std::vector<std::string>& modelVarToCalcList,
                   std::vector<std::string>& varRegrByList)
{
   std::string modelVarToCalcString;
   std::string varRegrByString;
   std::string filename;
   for (auto pathString : pathStrList)
   {
      modelVarToCalcString.clear();
      varRegrByString.clear();
      filename = ExtractFilename(pathString);
      ExtractModelVarForCalc(filename, modelVarToCalcString, varRegrByString);
      modelVarToCalcList.push_back(modelVarToCalcString);
      varRegrByList.push_back(varRegrByString);
   }
}


}  // namespace oops

#endif  // OOPS_GENERIC_STATVARIABLECHANGE_H_
