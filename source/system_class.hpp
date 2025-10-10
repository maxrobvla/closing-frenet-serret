#pragma once
#include <array>
#include <boost/math/interpolators/cubic_hermite.hpp>
#include <vector>

constexpr int dims = 3;
enum : unsigned int {
    POSITION = 0,
    TANGENT = dims,
    NORMAL = 2 * dims,
    BINORMAL = 3 * dims
};
typedef std::array<double, dims> real_vector;
typedef std::array<double, dims * dims> basis;
typedef std::array<double, (dims + 1) * dims> frenet_serret_frame;

struct harmonic_series_pure_cos {
    const int number_field_periods;
    std::vector<double> harmonic_coefficients;

    harmonic_series_pure_cos(int input_field_periods,
                             std::vector<double> input_coefficients);

    double operator()(double l);

#ifdef ENABLE_GRADIENT_CALCULATION
    double derivative_wrt_coefficient(double l, size_t coefficient_idx);
#endif
};

bool check_harmonic_series_for_zeros(int number_field_periods,
                                     std::vector<double> harmonic_coefficients);

struct general_curvature {
    const int number_field_periods;

    harmonic_series_pure_cos harmonic_part;

    const int order_zeros_field_maxima;
    const int order_zeros_field_minima;

    general_curvature(int input_field_periods,
                      std::vector<double> input_coefficients,
                      int order_zeros_field_maxima,
                      int order_zeros_field_minima);

    double zeros_func(double l);
    // helper function for operator()
    double half_period_function(double l);

    double operator()(double l);

#ifdef ENABLE_GRADIENT_CALCULATION
    double derivative_half_period_function(double l, size_t coefficient_idx);

    double derivative_wrt_coefficient(double l, size_t coefficient_idx);
#endif
};

typedef harmonic_series_pure_cos general_torsion;

#ifdef ENABLE_GRADIENT_CALCULATION
typedef boost::math::interpolators::cardinal_cubic_hermite_aos<
    std::vector<std::array<double, 2>>>
    interpolator;
#endif

class curve_system {
  public:
    size_t total_number_parameters; // total number of harmonic coefficients
    int number_field_periods;
    general_curvature curvature;
    general_torsion torsion;

    std::vector<frenet_serret_frame> curve;

    curve_system(int input_field_periods, std::vector<double> curvature_params,
                 std::vector<double> torsion_params,
                 int curvature_order_zeros_field_maxima,
                 int curvature_order_zeros_field_minima);

    void set_curve(bool prepare_derivative = false);

    void check_curve(bool prepare_derivative = false, bool verbose = false);

#ifdef ENABLE_GRADIENT_CALCULATION
    frenet_serret_frame interpolate(double l);

  private:
    std::vector<interpolator> interpolators;
#endif
};
