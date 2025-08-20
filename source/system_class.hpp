#pragma once
#include <array>
#include <boost/math/interpolators/cubic_hermite.hpp>
#include <vector>

constexpr int dims = 3;
typedef std::array<double, dims> real_vector;
typedef std::array<double, dims * dims> basis;
typedef std::array<double, (dims + 1) * dims> frenet_serret_frame;

struct harmonic_series_pure_cos {
    int number_field_periods;
    std::vector<double> harmonic_coefficients;

    harmonic_series_pure_cos(int input_field_periods,
                             std::vector<double> input_coefficients);

    double operator()(double l);

    double derivative_wrt_coefficient(double l, size_t coefficient_idx);
};

bool check_harmonic_series_for_zeros(int number_field_periods,
                                     std::vector<double> harmonic_coefficients);

struct general_curvature {
    int number_field_periods;

    harmonic_series_pure_cos harmonic_part;

    int order_zeros_field_maxima;
    int order_zeros_field_minima;

    general_curvature(int input_field_periods,
                      std::vector<double> input_coefficients,
                      int order_zeros_field_maxima,
                      int order_zeros_field_minima);

    double zeros_func(double l);
    // helper function for operator()
    double half_period_function(double l);

    double operator()(double l);

    double derivative_half_period_function(double l, size_t coefficient_idx);

    double derivative_wrt_coefficient(double l, size_t coefficient_idx);
};

typedef harmonic_series_pure_cos general_torsion;

typedef boost::math::interpolators::cardinal_cubic_hermite_aos<
    std::vector<std::array<double, 2>>>
    interpolator;

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

    void set_curve();

    void check_curve(bool verbose = false);

    frenet_serret_frame interpolate(double l);

  private:
    std::vector<interpolator> interpolators;
};
