#define sample_size (int)100
#define low_integration_bound (double)0
#define high_integration_bound (double)1
#define step_size                                                              \
    (double)((high_integration_bound - low_integration_bound) /                \
             (double)sample_size)
#include <array>

// we need half steps for integration
typedef std::array<double, sample_size> scalar_sample;
typedef std::array<double, 2 * sample_size> scalar_2x_sample;
typedef std::array<std::array<double, 3>, sample_size> vector_sample;
typedef std::array<double, 3> real_vector;
typedef std::array<std::array<double, 3>, 3> basis;
typedef std::array<basis, sample_size> basis_sample;

scalar_sample sample_function(double (*func)(double));
