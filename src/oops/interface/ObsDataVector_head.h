/*
 * (C) Copyright 2018  UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_OBSDATAVECTOR_HEAD_H_
#define OOPS_INTERFACE_OBSDATAVECTOR_HEAD_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace oops {
  template <typename OBS> class ObsSpace;
  template <typename OBS> class ObsVector;
  class ObsVariables;

// -----------------------------------------------------------------------------
/// \brief ObsDataVector is a vector templated on data type, in the observation space
/// \details
/// oops currently uses ObsDataVector<int> and ObsDataVector<float>);

template <typename OBS, typename DATATYPE>
class ObsDataVector : public util::Printable,
                      private util::ObjectCounter<ObsDataVector<OBS, DATATYPE> > {
  typedef typename OBS::template ObsDataVector<DATATYPE>  ObsDataVec_;

 public:
  static const std::string classname() {return "oops::ObsDataVector";}

  /// Constructor for specified ObsSpace \p os, with \p variables. If the group \p name is
  /// specified, the data is read from ObsSpace for specified variables and group.
  /// Otherwise ObsDataVector is allocated for specified variables and filled with zeros.
  ObsDataVector(const ObsSpace<OBS> & os, const ObsVariables & vars,
                         const std::string name = "");
  /// Copy constructor from \p other
  ObsDataVector(const ObsDataVector & other);
  /// Constructor from \p other ObsVector. ObsDataVector is created with variables from ObsVector
  /// and assigned ObsVector values. This is only well defined for numeric DATATYPE.
  explicit ObsDataVector(ObsVector<OBS> & other);
  /// Destructor (defined explicitly for timing and tracing)
  ~ObsDataVector();

  /// Accessor to the data
  ObsDataVec_ & obsdatavector() {return *data_;}
  /// const accessor to the data
  const ObsDataVec_ & obsdatavector() const {return *data_;}

  /// Accessor returning pointer to the data
  std::shared_ptr<ObsDataVec_> obsdatavectorptr() {return data_;}
  /// const accessor returning pointer to the data
  std::shared_ptr<const ObsDataVec_> obsdatavectorptr() const {return data_;}

  /// Assignment operator
  ObsDataVector & operator = (const ObsDataVector &);

  /// Zero out this ObsDataVector
  void zero();
  /// Mask values by reading another ObsDataVector \p qc that has the same variables and contains
  /// the masking information. Elements of *this* corresponding to non-zero elements of \p qc are
  /// set to missing.
  void mask(const ObsDataVector<OBS, int> & qc);
  /// Return the number of observations that aren't set to missing, across all MPI tasks.
  unsigned int nobs() const {return data_->nobs();}

  /// Fill ObsDataVector with data with group \p name from the associated ObsSpace
  void read(const std::string & name);
  /// Save this ObsDataVector as group \p name in the ObsSpace
  void save(const std::string & name) const;

 private:
  void print(std::ostream &) const;
  /// Pointer to the ObsDataVector implementation
  std::shared_ptr<ObsDataVec_> data_;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSDATAVECTOR_HEAD_H_
