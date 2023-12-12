/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_QGFORTRAN_H_
#define QG_MODEL_QGFORTRAN_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace qg {
  class LocationsQG;
  class ObsSpaceQG;

// Geometry key type
typedef int F90geom;
// Geometry iterator key type
typedef int F90iter;
// Model key type
typedef int F90model;
// Gom key type
typedef int F90gom;
// Fields key type
typedef int F90flds;
// Error covariance key type
typedef int F90error_covariance;
// Observation vector key type
typedef int F90ovec;
// Observation data base type
typedef int F90odb;
// Localization matrix
typedef int F90lclz;
// trajectory for wind speed
typedef int F90traj;

/// Interface to Fortran QG model
/*!
 * The core of the QG model is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {
// -----------------------------------------------------------------------------
//  Change of variable
// -----------------------------------------------------------------------------
  void qg_change_var_f90(const F90flds &, const oops::Variables &);
  void qg_change_var_tl_f90(const F90flds &, const oops::Variables &);
  void qg_change_var_ad_f90(const F90flds &, const oops::Variables &);

// -----------------------------------------------------------------------------
//  Error covariance
// -----------------------------------------------------------------------------
  void qg_error_covariance_setup_f90(F90error_covariance &, const eckit::Configuration &,
                                     const F90geom &);
  void qg_error_covariance_delete_f90(F90error_covariance &);
  void qg_error_covariance_mult_f90(const F90error_covariance &, const F90flds &, const F90flds &);
  void qg_error_covariance_inv_mult_f90(const F90error_covariance &, const F90flds &,
                                        const F90flds &);
  void qg_error_covariance_randomize_f90(const F90error_covariance &, const F90flds &);

// -----------------------------------------------------------------------------
//  Fields
// -----------------------------------------------------------------------------
  void qg_fields_create_f90(F90flds &, const F90geom &, const oops::Variables &, const bool &);
  void qg_fields_create_from_other_f90(F90flds &, const F90flds &, const F90geom &);
  void qg_fields_delete_f90(F90flds &);
  void qg_fields_zero_f90(const F90flds &);
  void qg_fields_ones_f90(const F90flds &);
  void qg_fields_dirac_f90(const F90flds &, const eckit::Configuration &);
  void qg_fields_random_f90(const F90flds &);
  void qg_fields_copy_f90(const F90flds &, const F90flds &);
  void qg_fields_self_add_f90(const F90flds &, const F90flds &);
  void qg_fields_self_sub_f90(const F90flds &, const F90flds &);
  void qg_fields_self_mul_f90(const F90flds &, const double &);
  void qg_fields_axpy_f90(const F90flds &, const double &, const F90flds &);
  void qg_fields_self_schur_f90(const F90flds &, const F90flds &);
  void qg_fields_dot_prod_f90(const F90flds &, const F90flds &, double &);
  void qg_fields_add_incr_f90(const F90flds &, const F90flds &);
  void qg_fields_diff_incr_f90(const F90flds &, const F90flds &, const F90flds &);
  void qg_fields_change_resol_f90(const F90flds &, const F90flds &);
  void qg_fields_change_resol_ad_f90(const F90flds &, const F90flds &);
  void qg_fields_read_file_f90(const F90flds &, const eckit::Configuration &,
                               util::DateTime &);
  void qg_fields_write_file_f90(const F90flds &, const eckit::Configuration &,
                                const util::DateTime &);
  void qg_fields_analytic_init_f90(const F90flds &, const eckit::Configuration &,
                                   util::DateTime &);
  void qg_fields_gpnorm_f90(const F90flds &, int[], double[], double[], double[]);
  void qg_fields_rms_f90(const F90flds &, double &);
  void qg_fields_sizes_f90(const F90flds &, int &, int &, int &);
  void qg_fields_lbc_f90(const F90flds &, int &);
  void qg_fields_to_fieldset_f90(const F90flds &,
                                 atlas::field::FieldSetImpl *);
  void qg_fields_from_fieldset_f90(const F90flds &,
                                   const atlas::field::FieldSetImpl *);
  void qg_fields_getvals_f90(const F90flds &, const oops::Variables &,
                             const int &, const double &, const int &, double &);
  void qg_fields_getvalsad_f90(const F90flds &, const oops::Variables &,
                               const int &, const double &, const int &, const double &);
  void qg_fields_getpoint_f90(const F90flds&, const F90iter&, const int &, double &);
  void qg_fields_setpoint_f90(const F90flds&, const F90iter&, const int &, const double &);
  void qg_fields_serialize_f90(const F90flds &, const int &, double[]);
  void qg_fields_deserialize_f90(const F90flds &, const int &, const double[], int &);

// -----------------------------------------------------------------------------
//  Geometry
// -----------------------------------------------------------------------------
  void qg_geom_setup_f90(F90geom &, const eckit::Configuration &);
  void qg_geom_set_lonlat_f90(const F90geom &, atlas::field::FieldSetImpl *);
  void qg_geom_set_functionspace_pointer_f90(const F90geom &,
                                                   atlas::functionspace::FunctionSpaceImpl *);
  void qg_geom_fill_geometry_fields_f90(const F90geom &, atlas::field::FieldSetImpl *);
  void qg_geom_clone_f90(F90geom &, const F90geom &);
  void qg_geom_info_f90(const F90geom &, int &, int &, int &, double &, double &);
  void qg_geom_delete_f90(F90geom &);
  void qg_geom_dimensions_f90(double &, double &, double &, double &, double &);

// -----------------------------------------------------------------------------
//  Geometry iterator
// -----------------------------------------------------------------------------
  void qg_geom_iter_setup_f90(F90iter &, const F90geom &, const int &);
  void qg_geom_iter_clone_f90(F90iter &, const F90iter &);
  void qg_geom_iter_delete_f90(F90iter &);
  void qg_geom_iter_equals_f90(const F90iter &, const F90iter&, int &);
  void qg_geom_iter_current_f90(const F90iter &, double &, double &);
  void qg_geom_iter_next_f90(const F90iter &);

// -----------------------------------------------------------------------------
//  Interpolation for getvalues via generic structures
// -----------------------------------------------------------------------------
  void qg_getvalues_interp_f90(const F90geom &, const atlas::field::FieldSetImpl *,
                               const oops::Variables &, const int &, const double &,
                               const int &, double &);
  void qg_getvalues_interp_ad_f90(const F90geom &, atlas::field::FieldSetImpl *,
                                  const oops::Variables &, const int &, const double &,
                                  const int &, const double &);

// -----------------------------------------------------------------------------
//  Local Values (GOM)
// -----------------------------------------------------------------------------
  void qg_gom_setup_f90(F90gom &, const int &, const oops::Variables &, const int &);
  void qg_gom_create_f90(F90gom &);
  void qg_gom_delete_f90(F90gom &);
  void qg_gom_fill_f90(const F90gom &, const int &, const char *, const int &, const int *,
                       const int &, const double *);
  void qg_gom_fillad_f90(const F90gom &, const int &, const char *, const int &, const int *,
                         const int &, double *);
  void qg_gom_copy_f90(const F90gom &, const F90gom &);
  void qg_gom_zero_f90(const F90gom &);
  void qg_gom_abs_f90(const F90gom &);
  void qg_gom_random_f90(const F90gom &);
  void qg_gom_mult_f90(const F90gom &, const double &);
  void qg_gom_add_f90(const F90gom &, const F90gom &);
  void qg_gom_diff_f90(const F90gom &, const F90gom &);
  void qg_gom_schurmult_f90(const F90gom &, const F90gom &);
  void qg_gom_divide_f90(const F90gom &, const F90gom &);
  void qg_gom_rms_f90(const F90gom &, double &);
  void qg_gom_dotprod_f90(const F90gom &, const F90gom &, double &);
  void qg_gom_stats_f90(const F90gom &, int &, double &, double &, double &);
  void qg_gom_maxloc_f90(const F90gom &, double &, int &, const oops::Variables &);
  void qg_gom_read_file_f90(const F90gom &, const oops::Variables &, const eckit::Configuration &);
  void qg_gom_write_file_f90(const F90gom &, const eckit::Configuration &);
  void qg_gom_analytic_init_f90(const F90gom &, const LocationsQG &,
                                const eckit::Configuration &);

// -----------------------------------------------------------------------------
//  Model
// -----------------------------------------------------------------------------
  void qg_model_setup_f90(F90model &, const eckit::Configuration &);
  void qg_model_delete_f90(F90model &);
  void qg_model_propagate_f90(const F90model &, const F90flds &);
  void qg_model_propagate_tl_f90(const F90model &, const F90flds &, const F90flds &);
  void qg_model_propagate_ad_f90(const F90model &, const F90flds &, const F90flds &);

// -----------------------------------------------------------------------------
//  Observation Handler
// -----------------------------------------------------------------------------
  void qg_obsdb_setup_f90(F90odb &, const eckit::Configuration &,
                          const util::DateTime &, const util::DateTime &);
  void qg_obsdb_delete_f90(F90odb &);
  void qg_obsdb_save_f90(const F90odb &);
  void qg_obsdb_get_f90(const F90odb &, const int &, const char *,
                        const int &, const char *, const F90ovec &);
  void qg_obsdb_put_f90(const F90odb &, const int &, const char *,
                        const int &, const char *, const F90ovec &);
  void qg_obsdb_locations_f90(const F90odb &, const int &, const char *,
                              atlas::field::FieldSetImpl *, std::vector<util::DateTime> &);
  void qg_obsdb_generate_f90(const F90odb &, const int &, const char *,
                             const eckit::Configuration &, const util::DateTime &,
                             const util::Duration &, const int &, int &);
  void qg_obsdb_nobs_f90(const F90odb &, const int &, const char *, int &);

// -----------------------------------------------------------------------------
//  Observation vector
// -----------------------------------------------------------------------------
  void qg_obsvec_setup_f90(F90ovec &, const int &, const int &);
  void qg_obsvec_clone_f90(F90ovec &, const F90ovec &);
  void qg_obsvec_delete_f90(F90ovec &);
  void qg_obsvec_copy_f90(const F90ovec &, const F90ovec &);
  void qg_obsvec_zero_f90(const F90ovec &);
  void qg_obsvec_settomissing_ith_f90(const F90ovec &, const int &);
  void qg_obsvec_ones_f90(const F90ovec &);
  /// set ObsVector (with key \p obsvector_key) values to missing values where
  /// mask ObsVector (with key \p mask_key) values are set to 1
  void qg_obsvec_mask_f90(const F90ovec & obsvector_key, const F90ovec & mask_key);
  /// set ObsVector (with key \p obsvector_key) values to missing values where
  /// mask ObsVector (with key \p mask_key) values are set to missing value
  void qg_obsvec_mask_with_missing_f90(const F90ovec & obsvector_key,
                                       const F90ovec & mask_key);
  void qg_obsvec_mul_scal_f90(const F90ovec &, const double &);
  void qg_obsvec_add_f90(const F90ovec &, const F90ovec &);
  void qg_obsvec_sub_f90(const F90ovec &, const F90ovec &);
  void qg_obsvec_mul_f90(const F90ovec &, const F90ovec &);
  void qg_obsvec_div_f90(const F90ovec &, const F90ovec &);
  void qg_obsvec_axpy_f90(const F90ovec &, const double &, const F90ovec &);
  void qg_obsvec_invert_f90(const F90ovec &);
  void qg_obsvec_random_f90(const ObsSpaceQG &, const F90ovec &);
  void qg_obsvec_dotprod_f90(const F90ovec &, const F90ovec &, double &);
  void qg_obsvec_stats_f90(const F90ovec &, double &, double &, double &);
  void qg_obsvec_nobs_f90(const F90ovec &, int &);
  void qg_obsvec_size_f90(const F90ovec &, int &);
  /// fill \p data (size \p nobs) with all non-masked out (non-missing) values
  void qg_obsvec_get_withmask_f90(const F90ovec &, const F90ovec & mask_key,
                                  double * data, const int & nobs);
  void qg_obsvec_nobs_withmask_f90(const F90ovec &, const F90ovec & mask_key, int &);


// -----------------------------------------------------------------------------
//  Streamfunction observations
// -----------------------------------------------------------------------------
  void qg_stream_equiv_f90(const F90odb &, const F90gom &, const F90ovec &, const double &);
  void qg_stream_equiv_tl_f90(const F90odb &, const F90gom &, const F90ovec &, const double &);
  void qg_stream_equiv_ad_f90(const F90odb &, const F90gom &, const F90ovec &, double &);

// -----------------------------------------------------------------------------
//  Wind observations
// -----------------------------------------------------------------------------
  void qg_wind_equiv_f90(const F90odb &, const F90gom &, const F90ovec &, const double &);
  void qg_wind_equiv_tl_f90(const F90odb &, const F90gom &, const F90ovec &, const double &);
  void qg_wind_equiv_ad_f90(const F90odb &, const F90gom &, const F90ovec &, double &);

// -----------------------------------------------------------------------------
//  Wind speed observations
// -----------------------------------------------------------------------------
  void qg_wspeed_equiv_f90(const F90odb &, const F90gom &, const F90ovec &, const double &);
  void qg_wspeed_equiv_tl_f90(const F90odb &, const F90gom &, const F90ovec &, const F90traj &,
                              const double &);
  void qg_wspeed_equiv_ad_f90(const F90odb &, const F90gom &, const F90ovec &, const F90traj &,
                              double &);
  void qg_wspeed_alloctraj_f90(const int &, const F90traj &);
  void qg_wspeed_settraj_f90(const F90odb &, const F90gom &, const F90traj &);
}

}  // namespace qg
#endif  // QG_MODEL_QGFORTRAN_H_
