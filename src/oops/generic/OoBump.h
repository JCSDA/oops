/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_OOBUMP_H_
#define OOPS_GENERIC_OOBUMP_H_

#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"
#include "oops/generic/oobump_f.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/util/Logger.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// OoBump C++ interface

class OoBump : private boost::noncopyable {
 public:
  OoBump(const UnstructuredGrid &, const eckit::LocalConfiguration,
         const eckit::mpi::Comm &,
         const int &, const int &, const int &, const int &);
  explicit OoBump(OoBump &);
  ~OoBump();

  int getKey() const {return keyOoBump_;}
  void deassociate() {keyOoBump_ = 0;}

  int getColocated() const;
  int getNts() const;
  int getCvSize() const;

  void addMember(const UnstructuredGrid &, const int &) const;
  void addPseudoMember(const UnstructuredGrid &, const int &) const;
  void removeMember(UnstructuredGrid &, const int &) const;
  void removePseudoMember(UnstructuredGrid &, const int &) const;
  void runDrivers() const;
  void multiplyVbal(UnstructuredGrid &) const;
  void multiplyVbalInv(UnstructuredGrid &) const;
  void multiplyVbalAd(UnstructuredGrid &) const;
  void multiplyVbalInvAd(UnstructuredGrid &) const;
  void multiplyNicas(UnstructuredGrid &) const;
  void randomizeNicas(UnstructuredGrid &) const;
  void getParam(const std::string &, UnstructuredGrid &) const;
  void setParam(const std::string &, const UnstructuredGrid &) const;

 private:
  int keyOoBump_;
};

// -----------------------------------------------------------------------------
OoBump::OoBump(const UnstructuredGrid & ug, const eckit::LocalConfiguration conf,
               const eckit::mpi::Comm & comm,
               const int & ens1_ne, const int & ens1_nsub,
               const int & ens2_ne, const int & ens2_nsub) : keyOoBump_(0) {
  const eckit::Configuration * fconf = &conf;
  oobump_create_f90(keyOoBump_, ug.toFortran(), &fconf, ens1_ne, ens1_nsub,
                    ens2_ne, ens2_nsub, comm.name().size(), comm.name().c_str());
}
// -----------------------------------------------------------------------------
OoBump::OoBump(OoBump & other) : keyOoBump_(other.getKey()) {
  other.deassociate();
}
// -----------------------------------------------------------------------------
OoBump::~OoBump() {
  if (keyOoBump_ > 0) oobump_delete_f90(keyOoBump_);
}
// -----------------------------------------------------------------------------

int OoBump::getColocated() const {
  int colocated;
  oobump_get_colocated_f90(keyOoBump_, colocated);
  return colocated;
}
// -----------------------------------------------------------------------------
int OoBump::getNts() const {
  int nts;
  oobump_get_nts_f90(keyOoBump_, nts);
  return nts;
}
// -----------------------------------------------------------------------------
int OoBump::getCvSize() const {
  int cv_size;
  oobump_get_cv_size_f90(keyOoBump_, cv_size);
  return cv_size;
}
// -----------------------------------------------------------------------------
void OoBump::addMember(const UnstructuredGrid & ug, const int & ie) const {
  oobump_add_member_f90(keyOoBump_, ug.toFortran(), ie+1, 1);
}
// -----------------------------------------------------------------------------
void OoBump::addPseudoMember(const UnstructuredGrid & ug, const int & ie) const {
  oobump_add_member_f90(keyOoBump_, ug.toFortran(), ie+1, 2);
}
// -----------------------------------------------------------------------------
void OoBump::removeMember(UnstructuredGrid & ug, const int & ie) const {
  oobump_remove_member_f90(keyOoBump_, ug.toFortran(), ie+1, 1);
}
// -----------------------------------------------------------------------------
void OoBump::removePseudoMember(UnstructuredGrid & ug, const int & ie) const {
  oobump_remove_member_f90(keyOoBump_, ug.toFortran(), ie+1, 2);
}
// -----------------------------------------------------------------------------
void OoBump::runDrivers() const {
  oobump_run_drivers_f90(keyOoBump_);
}
// -----------------------------------------------------------------------------
void OoBump::multiplyVbal(UnstructuredGrid & ug) const {
  oobump_multiply_vbal_f90(keyOoBump_, ug.toFortran());
}
// -----------------------------------------------------------------------------
void OoBump::multiplyVbalInv(UnstructuredGrid & ug) const {
  oobump_multiply_vbal_inv_f90(keyOoBump_, ug.toFortran());
}
// -----------------------------------------------------------------------------
void OoBump::multiplyVbalAd(UnstructuredGrid & ug) const {
  oobump_multiply_vbal_ad_f90(keyOoBump_, ug.toFortran());
}
// -----------------------------------------------------------------------------
void OoBump::multiplyVbalInvAd(UnstructuredGrid & ug) const {
  oobump_multiply_vbal_inv_ad_f90(keyOoBump_, ug.toFortran());
}
// -----------------------------------------------------------------------------
void OoBump::multiplyNicas(UnstructuredGrid & ug) const {
  oobump_multiply_nicas_f90(keyOoBump_, ug.toFortran());
}
// -----------------------------------------------------------------------------
void OoBump::randomizeNicas(UnstructuredGrid & ug) const {
  oobump_randomize_nicas_f90(keyOoBump_, ug.toFortran());
}
// -----------------------------------------------------------------------------
void OoBump::getParam(const std::string & param, UnstructuredGrid & ug) const {
  const int nstr = param.size();
  const char *cstr = param.c_str();
  oobump_get_param_f90(keyOoBump_, nstr, cstr, ug.toFortran());
}
// -----------------------------------------------------------------------------
void OoBump::setParam(const std::string & param, const UnstructuredGrid & ug) const {
  const int nstr = param.size();
  const char *cstr = param.c_str();
  oobump_set_param_f90(keyOoBump_, nstr, cstr, ug.toFortran());
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_OOBUMP_H_
