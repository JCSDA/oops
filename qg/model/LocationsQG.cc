#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "model/QgFortran.h"
#include "model/LocationsQG.h"

namespace qg {

  // -------------------------------------------------------------------------
  /*! QG Locations Constructor with Configuration
   * 
   * \details This constructor can be used to generate user-specified
   * and/or random locations for use with testStateInterpolation()
   *
   * To generate random locations, the relevant parameters specified in 
   * **State.Locations** section of the config file are:
   * 
   * * **lats** user-specified latitudes (degrees)
   * * **lons** user-specified longitudes (degrees)
   * * **heights** user-specified heights (normalized, from 0-1)
   * * **Nrandom** number of random locations desired
   * * **random_seed** (optional) random seed for reproducibility of results
   * 
   * \date April, 2018 Created (M. Miesch, JCSDA)
   *
   * \warning latitudes and longitudes are converted to normalized locations 
   * (between 0-1) based on on an assumed latitudinal extent of 180 degrees
   * and an assumed longitudinal extent of 360 degrees
   *
   * \sa qg::c_qg_loc_test() test::testStateInterpolation() 
   *
   */
 
  LocationsQG::LocationsQG(const eckit::Configuration & config) {
    const eckit::Configuration * conf = &config;

    qg_loc_create_f90(keyLoc_);

    if (config.has("lats") || config.has("Nrandom")) {
	std::vector<double> lats = config.getDoubleVector("lats");
	std::vector<double> lons = config.getDoubleVector("lons");
	std::vector<double> zin  = config.getDoubleVector("heights");

	ASSERT(lats.size() == lons.size());
	const unsigned int nloc = lats.size();

	// Default to level 1 unless otherwise specified
	std::vector<double> height(nloc,0.25);
	
	if (zin.size() > 0) {
	  for (unsigned int i=0; i < zin.size(); ++i) {
	    if (i >= nloc)
	      break;
	    height[i] = zin[i];
	  }
	}
	
	qg_loc_test_f90(keyLoc_,&conf,nloc,&lats[0],&lons[0],&height[0]);
      }
  }

  // -------------------------------------------------------------------------
}
