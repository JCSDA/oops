/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_QGFORTRAN_H_
#define QG_MODEL_QGFORTRAN_H_

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace qg {

// Geometry key type
typedef int F90geom;
// Model key type
typedef int F90model;
// Variables key type
typedef int F90vars;
// Locations key type
typedef int F90locs;
// Goms key type
typedef int F90goms;
// Fields key type
typedef int F90flds;
// Trajectory key type
typedef int F90traj;
// Background error covariance key type
typedef int F90bmat;
// Observation vector key type
typedef int F90ovec;
// Obs operator key type
typedef int F90hop;
// Observation data base type
typedef int F90odb;
// Localization matrix
typedef int F90lclz;

/// Interface to Fortran QG model
/*!
 * The core of the QG model is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {
// -----------------------------------------------------------------------------
//  Geometry
// -----------------------------------------------------------------------------
  void qg_geo_setup_f90(F90geom &, const eckit::Configuration * const *);
  void qg_geo_clone_f90(const F90geom &, F90geom &);
  void qg_geo_info_f90(const F90geom &, int &, int &);
  void qg_geo_delete_f90(F90geom &);

// -----------------------------------------------------------------------------
//  Model
// -----------------------------------------------------------------------------
  void qg_setup_f90(const eckit::Configuration * const *, const F90geom &, F90model &);
  void qg_delete_f90(F90model &);

  void qg_prepare_integration_f90(const F90model &, const F90flds &);
  void qg_prepare_integration_tl_f90(const F90model &, const F90flds &);
  void qg_prepare_integration_ad_f90(const F90model &, const F90flds &);

  void qg_propagate_f90(const F90model &, const F90flds &);
  void qg_prop_traj_f90(const F90model &, const F90flds &, F90traj &);
  void qg_propagate_tl_f90(const F90model &, const F90flds &, const F90traj &);
  void qg_propagate_ad_f90(const F90model &, const F90flds &, const F90traj &);

  void qg_wipe_traj_f90(F90traj &);
  void qg_traj_minmaxrms_f90(const F90traj &, double &);

// -----------------------------------------------------------------------------
//  Fields
// -----------------------------------------------------------------------------
  void qg_field_create_f90(F90flds &, const F90geom &, const F90vars *);
  void qg_field_delete_f90(F90flds &);

  void qg_field_copy_f90(const F90flds &, const F90flds &);
  void qg_field_zero_f90(const F90flds &);
  void qg_field_self_add_f90(const F90flds &, const F90flds &);
  void qg_field_self_sub_f90(const F90flds &, const F90flds &);
  void qg_field_self_mul_f90(const F90flds &, const double &);
  void qg_field_axpy_f90(const F90flds &, const double &, const F90flds &);
  void qg_field_dot_prod_f90(const F90flds &, const F90flds &, double &);
  void qg_field_self_schur_f90(const F90flds &, const F90flds &);
  void qg_field_random_f90(const F90flds &);

  void qg_field_add_incr_f90(const F90flds &, const F90flds &);
  void qg_field_diff_incr_f90(const F90flds &, const F90flds &, const F90flds &);

  void qg_field_change_resol_f90(const F90flds &, const F90flds &);

  void qg_field_read_file_f90(const F90flds &, const eckit::Configuration * const *,
                              util::DateTime * const *);
  void qg_field_analytic_init_f90(const F90flds &, const F90geom &,
                                  const eckit::Configuration * const *,
                                  util::DateTime * const *);
  void qg_field_write_file_f90(const F90flds &, const eckit::Configuration * const *,
                               const util::DateTime * const *);

  void qg_field_interp_f90(const F90flds &, const F90locs &, const F90vars *, const F90goms &);
  void qg_field_interp_tl_f90(const F90flds &, const F90locs &, const F90vars *, const F90goms &);
  void qg_field_interp_ad_f90(const F90flds &, const F90locs &, const F90vars *, const F90goms &);

  void qg_field_gpnorm_f90(const F90flds &, const int &, double &);
  void qg_field_sizes_f90(const F90flds &, int &, int &, int &, int &);
  void qg_field_rms_f90(const F90flds &, double &);

// -----------------------------------------------------------------------------
//  Background error
// -----------------------------------------------------------------------------
  void qg_b_setup_f90(F90bmat &, const eckit::Configuration * const *, const F90geom &);
  void qg_b_delete_f90(F90bmat &);

  void qg_b_linearize_f90(const F90bmat &, const eckit::Configuration * const *);

  void qg_b_mult_f90(const F90bmat &, const F90flds &, const F90flds &);
  void qg_b_invmult_f90(const F90bmat &, const F90flds &, const F90flds &);

  void qg_b_randomize_f90(const F90bmat &, const F90flds &);

// -----------------------------------------------------------------------------
//  Localization matrix
// -----------------------------------------------------------------------------
  void qg_localization_setup_f90(F90lclz &, const eckit::Configuration * const *,
                                 const F90geom &);
  void qg_localization_delete_f90(F90lclz &);
  void qg_localization_mult_f90(const F90lclz &, const F90flds &);

// -----------------------------------------------------------------------------
//  Locations
// -----------------------------------------------------------------------------
  void qg_loc_create_f90(F90locs &);
  void qg_loc_test_f90(const F90locs &, const eckit::Configuration * const *,
                       const int &, const double *, const double *, const double*);
  void qg_loc_delete_f90(F90locs &);
  void qg_loc_nobs_f90(const F90locs &, int &);
  void qg_loc_element_f90(const F90locs &, const int &, double *);

// -----------------------------------------------------------------------------
//  Local Values (GOM)
// -----------------------------------------------------------------------------
  void qg_gom_setup_f90(F90goms &, const F90locs &, const F90vars *);
  void qg_gom_create_f90(F90goms &);
  void qg_gom_delete_f90(F90goms &);
  void qg_gom_abs_f90(const F90goms &);
  void qg_gom_rms_f90(const F90goms &, double &);
  void qg_gom_zero_f90(const F90goms &);
  void qg_gom_assign_f90(const F90goms &, const F90goms &);
  void qg_gom_random_f90(const F90goms &);
  void qg_gom_mult_f90(const F90goms &, const double &);
  void qg_gom_add_f90(const F90goms &, const F90goms &);
  void qg_gom_diff_f90(const F90goms &, const F90goms &);
  void qg_gom_normalize_f90(const F90goms &, const F90goms &);
  void qg_gom_dotprod_f90(const F90goms &, const F90goms &, double &);
  void qg_gom_minmaxavg_f90(const F90goms &, int &, double &, double &, double &);
  void qg_gom_maxloc_f90(const F90goms &, double &, int &, int &);
  void qg_gom_read_file_f90(const F90goms &, const eckit::Configuration * const *);
  void qg_gom_analytic_init_f90(const F90goms &, const F90locs &,
                                const eckit::Configuration * const *);
  void qg_gom_write_file_f90(const F90goms &, const eckit::Configuration * const *);

// -----------------------------------------------------------------------------
//  Streamfunction observations
// -----------------------------------------------------------------------------
  void qg_stream_setup_f90(F90hop &, const eckit::Configuration * const *);
  void qg_stream_delete_f90(F90hop &);

  void qg_stream_equiv_f90(const F90goms &, const F90ovec &, const double &);
  void qg_stream_equiv_tl_f90(const F90goms &, const F90ovec &, const double &);
  void qg_stream_equiv_ad_f90(const F90goms &, const F90ovec &, double &);

// -----------------------------------------------------------------------------
//  Wind observations
// -----------------------------------------------------------------------------
  void qg_wind_setup_f90(F90hop &, const eckit::Configuration * const *);
  void qg_wind_delete_f90(F90hop &);

  void qg_wind_equiv_f90(const F90goms &, F90ovec &, const double &);
  void qg_wind_equiv_tl_f90(const F90goms &, const F90ovec &, const double &);
  void qg_wind_equiv_ad_f90(const F90goms &, const F90ovec &, double &);

// -----------------------------------------------------------------------------
//  Wind speed observations
// -----------------------------------------------------------------------------
  void qg_wspeed_setup_f90(F90hop &, const eckit::Configuration * const *);
  void qg_wspeed_delete_f90(F90hop &);

  void qg_wspeed_eqv_f90(const F90goms &, const F90ovec &, const double &);
  void qg_wspeed_equiv_tl_f90(const F90goms &, const F90ovec &, const F90goms &, const double &);
  void qg_wspeed_equiv_ad_f90(const F90goms &, const F90ovec &, const F90goms &, double &);

  void qg_wspeed_gettraj_f90(const F90hop &, const int &, const int *, F90goms &);
  void qg_wspeed_settraj_f90(const F90goms &, const F90goms &);

// -----------------------------------------------------------------------------
//  Observation Vectors
// -----------------------------------------------------------------------------
  void qg_obsvec_setup_f90(F90ovec &, const int &, const int &);
  void qg_obsvec_clone_f90(const F90ovec &, F90ovec &);
  void qg_obsvec_delete_f90(F90ovec &);

  void qg_obsvec_assign_f90(const F90ovec &, const F90ovec &);
  void qg_obsvec_zero_f90(const F90ovec &);
  void qg_obsvec_mul_scal_f90(const F90ovec &, const double &);
  void qg_obsvec_add_f90(const F90ovec &, const F90ovec &);
  void qg_obsvec_sub_f90(const F90ovec &, const F90ovec &);
  void qg_obsvec_mul_f90(const F90ovec &, const F90ovec &);
  void qg_obsvec_div_f90(const F90ovec &, const F90ovec &);
  void qg_obsvec_axpy_f90(const F90ovec &, const double &, const F90ovec &);
  void qg_obsvec_invert_f90(const F90ovec &);
  void qg_obsvec_random_f90(const F90ovec &);
  void qg_obsvec_dotprod_f90(const F90ovec &, const F90ovec &, double &);
  void qg_obsvec_minmaxavg_f90(const F90ovec &, double &, double &, double &);
  void qg_obsvec_nobs_f90(const F90ovec &, int &);

// -----------------------------------------------------------------------------
//  Observation Handler
// -----------------------------------------------------------------------------
  void qg_obsdb_setup_f90(F90odb &, const eckit::Configuration * const *);
  void qg_obsdb_delete_f90(F90odb &);
  void qg_obsdb_get_f90(const F90odb &, const int &, const char *,
                        const int &, const char *, const F90ovec &);
  void qg_obsdb_put_f90(const F90odb &, const int &, const char *,
                        const int &, const char *, const F90ovec &);
  void qg_obsdb_locations_f90(const F90odb &, const int &, const char *,
                              const util::DateTime * const *, const util::DateTime * const *,
                              F90locs &);
  void qg_obsdb_generate_f90(const F90odb &, const int &, const char *,
                             const eckit::Configuration * const *, const util::DateTime * const *,
                             const util::Duration * const *, const int &, int &);
  void qg_obsdb_nobs_f90(const F90odb &, const int &, const char *, int &);
  void qg_obsoper_inputs_f90(const F90hop &, F90vars *);

}

// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_QGFORTRAN_H_
