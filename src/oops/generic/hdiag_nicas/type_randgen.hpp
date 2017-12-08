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

        // Random integer generator
        void rand_integer(int binf, int bsup, int *ir) const;

        // Random real generator
        void rand_real(double binf, double bsup, double *rr) const;

        // Sampling initialization
        void initialize_sampling(int n, double lon[], double lat[], int mask[], double rh[], int ntry, int nrep, int ns, int ihor[]) const;
    private:
#if __cplusplus > 199711L
        // Mersenne Twister 19937 generator
        std::mt19937 *gen;
#endif
};
