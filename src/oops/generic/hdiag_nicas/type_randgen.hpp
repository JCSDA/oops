//---------------------------------------------------------------------
/// Purpose: random numbers generator class header
/// Author : Benjamin Menetrier
/// Licensing: this code is distributed under the CeCILL-C license
/// Copyright © 2017 METEO-FRANCE
// ----------------------------------------------------------------------
#if __cplusplus > 199711L
#include <random>
#endif

class randGen {
    public:
        // Constructor
        randGen(int default_seed);

        // Destructor
        ~randGen();

        // Get version
        void get_version(int *version) const {*version=version_;};

        // Random integer generator
        void rand_integer(int binf, int bsup, int *ir);

        // Random real generator
        void rand_real(double binf, double bsup, double *rr);

        // Sampling initialization
        void initialize_sampling(int n, double lon[], double lat[], int mask[], double rh[], int ntry, int nrep, int ns, int ihor[]);
    private:
        // Version key
        int version_;

#if __cplusplus > 199711L
        // Mersenne Twister 19937 generator
        std::mt19937 *gen;
#endif

        // ran3 generator
        int seed;
        double ran3();
};
