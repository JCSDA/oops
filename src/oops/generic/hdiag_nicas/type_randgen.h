//---------------------------------------------------------------------
/// Purpose: random numbers generator class header
/// Author : Benjamin Menetrier
/// Licensing: this code is distributed under the CeCILL-C license
/// Copyright © 2017 METEO-FRANCE
// ----------------------------------------------------------------------
extern "C"
{
    class randGen;
    typedef randGen RANDGEN;

    // Constructor
    RANDGEN* create_randgen(unsigned long int default_seed);

    // Destructor
    void delete_randgen(RANDGEN* randgen);

    // Reseed generator
    void reseed_randgen(RANDGEN* randgen, unsigned long int seed);

    // Get version
    void get_version(const RANDGEN* randgen, int *version);

    // Random integer generator
    void rand_integer(RANDGEN* randgen, int binf, int bsup, int *ir);

    // Random real generator
    void rand_real(RANDGEN* randgen, double binf, double bsup, double *rr);

    // Sampling initialization
    void initialize_sampling(RANDGEN* randgen, int n, double lon[], double lat[], int mask[], double rh[], int ntry, int nrep, int ns, int ihor[]);
}
