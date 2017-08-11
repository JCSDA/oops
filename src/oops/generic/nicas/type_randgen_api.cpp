//---------------------------------------------------------------------
/// Purpose: random number generator class API
/// Author : Benjamin Menetrier
/// Licensing: this code is distributed under the CeCILL-C license
/// Copyright © 2017 METEO-FRANCE
// ----------------------------------------------------------------------
#include "type_randgen.h"
#include "type_randgen.hpp"

using namespace std;

// Constructor
RANDGEN* create_randgen(int default_seed){
    return new randGen(default_seed);
}

// Destructor
void delete_randgen(RANDGEN* randgen){
    delete randgen;
}

// Random integer generator
void rand_integer(const RANDGEN* randgen, int binf, int bsup, int *ir) {
    randgen->rand_integer(binf, bsup, ir);
}

// Sampling initialization
void initialize_sampling(const RANDGEN* randgen, int n, double lon[], double lat[], int mask[], double L[], int ntry, int nrep, int ns, int nfor, int ifor[], int ihor[]) {
    randgen->initialize_sampling( n, lon, lat, mask, L, ntry, nrep, ns, nfor, ifor, ihor);
}
