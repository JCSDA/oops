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
    RANDGEN* create_randgen(int default_seed);

    // Destructor
    void delete_randgen(RANDGEN* randgen);

    // Random integer generator
    void rand_integer(const RANDGEN* randgen, int binf, int bsup, int *ir);

    // Random real generator
    void rand_real(const RANDGEN* randgen, double binf, double bsup, double *rr);

    // Sampling initialization
    void initialize_sampling(const RANDGEN* randgen, int n, double lon[], double lat[], int mask[], double L[], int ntry, int nrep, int ns, int nfor, int ifor[], int ihor[]);
}
