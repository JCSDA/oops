/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_OOBUMP_H_
#define OOPS_GENERIC_OOBUMP_H_

#include <algorithm>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

#include "oops/assimilation/Increment4D.h"
#include "oops/base/Variables.h"
#include "oops/generic/oobump_f.h"
#if !ATLASIFIED
#include "oops/generic/UnstructuredGrid.h"
#endif
#include "oops/interface/Increment.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// OoBump C++ interface

template<typename MODEL> class OoBump : public boost::noncopyable {
  typedef Geometry<MODEL>    Geometry_;
  typedef Increment<MODEL>   Increment_;
  typedef Increment4D<MODEL> Increment4D_;

 public:
  OoBump(const Geometry_ &, const Variables &, const std::vector<util::DateTime> &,
         const eckit::LocalConfiguration, const int &, const int &, const int &, const int &);
  explicit OoBump(OoBump &);
  ~OoBump();

  // C++ interfaces
  size_t getSize() {return keyOoBump_.size();}
  int getKey(int igrid) const {return keyOoBump_[igrid];}
  void clearKey() {keyOoBump_.clear();}

  // Fortran interfaces
  int getCvSize() const;
  void addMember(Increment4D_ &, const int &) const;
  void addPseudoMember(Increment4D_ &, const int &) const;
  void removeMember(Increment4D_ &, const int &) const;
  void removePseudoMember(Increment4D_ &, const int &) const;
  void runDrivers() const;
  void multiplyVbal(const Increment_ &, Increment_ &) const;
  void multiplyVbalInv(const Increment_ &, Increment_ &) const;
  void multiplyVbalAd(const Increment_ &, Increment_ &) const;
  void multiplyVbalInvAd(const Increment_ &, Increment_ &) const;
  void multiplyNicas(Increment_ &) const;
  void multiplyNicas(const Increment_ &, Increment_ &) const;
  void multiplyNicas(Increment4D_ &) const;
  void multiplyNicas(const Increment4D_ &, Increment4D_ &) const;
  void randomizeNicas(Increment_ &) const;
  void randomizeNicas(Increment4D_ &) const;
  void getParam(const std::string &, Increment4D_ &) const;
  void setParam(const std::string &, const Increment4D_ &) const;

 private:
  std::vector<int> keyOoBump_;
#if !ATLASIFIED
  std::unique_ptr<UnstructuredGrid> ug_;
#endif
};

// -----------------------------------------------------------------------------
template<typename MODEL>
OoBump<MODEL>::OoBump(const Geometry_ & resol,
                      const Variables & vars,
                      const std::vector<util::DateTime> & timeslots,
                      const eckit::LocalConfiguration conf,
                      const int & ens1_ne, const int & ens1_nsub,
                      const int & ens2_ne, const int & ens2_nsub) : keyOoBump_() {
  // Get configuration pointer
  const eckit::Configuration * fconf = &conf;

  // Grids
  std::vector<eckit::LocalConfiguration> grids;

#if ATLASIFIED
  // Get the grids configuration from input configuration and complete it

  if (conf.has("grids")) {
    // Get grids from input configuration
    conf.get("grids", grids);
    ASSERT(grids.size() > 0);
  } else {
    // Create one empty configuration
    eckit::LocalConfiguration emptyConf;
    grids.push_back(emptyConf);
  }

  // Loop over grids
  for (unsigned int jgrid = 0; jgrid < grids.size(); ++jgrid) {
    // Add grid index
    int jgrid_int(jgrid);
    grids[jgrid].set("grid_index", jgrid_int);

    // Add input variables to the grid configuration
    std::vector<std::string> vars_str;
    if (grids[jgrid].has("variables")) {
      grids[jgrid].get("variables", vars_str);
    } else {
      for (unsigned int jvar = 0; jvar < vars.size(); ++jvar) {
        vars_str.push_back(vars[jvar]);
      }
      grids[jgrid].set("variables", vars_str);
    }

    // Add input timeslots to the grid configuration
    std::vector<std::string> timeslots_str;
    if (grids[jgrid].has("timeslots")) {
      grids[jgrid].get("timeslots", timeslots_str);
    } else {
      for (unsigned int jts = 0; jts < timeslots.size(); ++jts) {
        timeslots_str.push_back(timeslots[jts].toString());
      }
      grids[jgrid].set("timeslots", timeslots_str);
    }

    // Get the required number of levels add it to the grid configuration
    Increment_ dx(resol, vars, timeslots[0]);
    std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
    dx.setAtlas(atlasFieldSet.get());
    int nl = 0;
    for (unsigned int jvar = 0; jvar < vars_str.size(); ++jvar) {
      std::string fieldname = vars_str[jvar] + '_' + timeslots_str[0];
      atlas::Field atlasField = atlasFieldSet->field(fieldname);
      nl = std::max(nl, atlasField.levels());
    }
    grids[jgrid].set("nl", nl);

    // Add level index for 2D fields (first or last, first by default)
    if (!grids[jgrid].has("lev2d")) {
      grids[jgrid].set("lev2d", "first");
    }
  }
#else
  // Get the grids configuration from the unstructured grid configuration
  Increment4D_ dx(resol, vars, timeslots);
  int colocated = 1;
  if (conf.has("colocated")) {
    colocated = conf.getInt("colocated");
  }
  ug_.reset(new UnstructuredGrid(colocated, timeslots.size()));
  dx.ug_coord(*ug_.get());
  ug_->defineGeometry();
  ug_->defineGrids(grids);
#endif

// Check grids number
  ASSERT(grids.size() > 0);

  for (unsigned int jgrid = 0; jgrid < grids.size(); ++jgrid) {
    // Print configuration for this grid
    Log::info() << "Grid " << jgrid << ": " << grids[jgrid] << std::endl;

    // Get grid configuration pointer
    const eckit::Configuration * fgrid = &grids[jgrid];

    // Create OoBump instance
    int keyOoBump = 0;
#if ATLASIFIED
    oobump_create_f90(keyOoBump, &resol.getComm(),
                      resol.atlasFunctionSpace()->get(),
                      resol.atlasFieldSet()->get(),
                      &fconf, &fgrid,
                      ens1_ne, ens1_nsub, ens2_ne, ens2_nsub);
#else
    oobump_create_f90(keyOoBump, &resol.getComm(),
                      ug_->atlasFunctionSpace()->get(),
                      ug_->atlasFieldSet()->get(),
                      &fconf, &fgrid,
                      ens1_ne, ens1_nsub, ens2_ne, ens2_nsub);
#endif
    keyOoBump_.push_back(keyOoBump);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
OoBump<MODEL>::OoBump(OoBump & other) : keyOoBump_() {
  for (unsigned int jgrid = 0; jgrid < other.getSize(); ++jgrid) {
    keyOoBump_.push_back(other.getKey(jgrid));
  }
  other.clearKey();
#if !ATLASIFIED
  ug_ = std::move(other.ug_);
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
OoBump<MODEL>::~OoBump() {
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    if (keyOoBump_[jgrid] > 0) oobump_delete_f90(keyOoBump_[jgrid]);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
int OoBump<MODEL>::getCvSize() const {
  int cv_size_tot = 0;
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    int cv_size = 0;
    oobump_get_cv_size_f90(keyOoBump_[jgrid], cv_size);
    cv_size_tot += cv_size;
  }
  return cv_size_tot;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::addMember(Increment4D_ & dx, const int & ie) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.toAtlas(atlasFieldSet.get());
#else
  dx.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    oobump_add_member_f90(keyOoBump_[jgrid], atlasFieldSet->get(), ie+1, 1);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::addPseudoMember(Increment4D_ & dx, const int & ie) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.toAtlas(atlasFieldSet.get());
#else
  dx.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    oobump_add_member_f90(keyOoBump_[jgrid], atlasFieldSet->get(), ie+1, 2);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::removeMember(Increment4D_ & dx, const int & ie) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.setAtlas(atlasFieldSet.get());
#else
  ug_->setAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    oobump_remove_member_f90(keyOoBump_[jgrid], atlasFieldSet->get(), ie+1, 1);
  }
#if ATLASIFIED
  dx.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dx.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::removePseudoMember(Increment4D_ & dx, const int & ie) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.setAtlas(atlasFieldSet.get());
#else
  ug_->setAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    oobump_remove_member_f90(keyOoBump_[jgrid], atlasFieldSet->get(), ie+1, 2);
  }
#if ATLASIFIED
  dx.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dx.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::runDrivers() const {
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    oobump_run_drivers_f90(keyOoBump_[jgrid]);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyVbal(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    oobump_multiply_vbal_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyVbalInv(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    oobump_multiply_vbal_inv_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyVbalAd(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    oobump_multiply_vbal_ad_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyVbalInvAd(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    oobump_multiply_vbal_inv_ad_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyNicas(Increment_ & dx) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.toAtlas(atlasFieldSet.get());
#else
  dx.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    oobump_multiply_nicas_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dx.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dx.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyNicas(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    oobump_multiply_nicas_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyNicas(Increment4D_ & dx) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.toAtlas(atlasFieldSet.get());
#else
  dx.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    oobump_multiply_nicas_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dx.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dx.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyNicas(const Increment4D_ & dxi, Increment4D_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    oobump_multiply_nicas_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::randomizeNicas(Increment_ & dx) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.setAtlas(atlasFieldSet.get());
#else
  ug_->setAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    oobump_randomize_nicas_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dx.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dx.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::randomizeNicas(Increment4D_ & dx) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.setAtlas(atlasFieldSet.get());
#else
  ug_->setAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    oobump_randomize_nicas_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dx.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dx.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::getParam(const std::string & param, Increment4D_ & dx) const {
  const int nstr = param.size();
  const char *cstr = param.c_str();
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.setAtlas(atlasFieldSet.get());
#else
  ug_->setAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    oobump_get_param_f90(keyOoBump_[jgrid], nstr, cstr, atlasFieldSet.get()->get());
  }
#if ATLASIFIED
  dx.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dx.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::setParam(const std::string & param, const Increment4D_ & dx) const {
  const int nstr = param.size();
  const char *cstr = param.c_str();
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.toAtlas(atlasFieldSet.get());
#else
  dx.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    oobump_set_param_f90(keyOoBump_[jgrid], nstr, cstr, atlasFieldSet->get());
  }
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_OOBUMP_H_
