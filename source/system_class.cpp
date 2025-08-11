#include "system_class.hpp"
#include "solver.hpp"
#include <cmath>
#include <vector>

harmonic_series_pure_cos::harmonic_series_pure_cos(
    int input_field_periods, std::vector<double> input_coefficients)
    : number_field_periods(input_field_periods),
      harmonic_coefficients(input_coefficients) {};

double harmonic_series_pure_cos::operator()(double l) {
    double result = harmonic_coefficients[0];
    for (size_t i = 1; i < harmonic_coefficients.size(); i++) {
        result += harmonic_coefficients[i] *
                  std::cos(i * double(number_field_periods) * l);
    }
    return result;
}

bool check_harmonic_series_for_zeros(harmonic_series_pure_cos series) {
    double half_period_length =
        std::numbers::pi / double(series.number_field_periods);
    double DX = half_period_length / double(sample_size);
    double l;

    std::array<double, sample_size> grid_values;
    for (int i = 0; i < sample_size; i++) {
        l = DX * double(i);
        grid_values[i] = series(l);
    }

    // check the sign of every component except the first and the last, because
    // those are 0 anyways
    // idea: if the sign changes throughout the half field period, there are
    // zeros inside the half field period
    bool only_zeros_on_boundary = true;
    bool positive = (grid_values[1] > 0);
    for (size_t i = 2; i < sample_size - 1; i++) {
        if ((positive and (grid_values[i] < 0)) or
            ((not positive) and (grid_values[i] > 0))) {
            only_zeros_on_boundary = false;
            break;
        }
    }
    return only_zeros_on_boundary;
}

bool check_harmonic_series_for_zeros(
    int number_field_periods, std::vector<double> harmonic_coefficients) {
    harmonic_series_pure_cos series(number_field_periods,
                                    harmonic_coefficients);
    return check_harmonic_series_for_zeros(series);
}

general_curvature::general_curvature(int input_field_periods,
                                     std::vector<double> input_coefficients,
                                     int input_order_zeros_field_maxima,
                                     int input_order_zeros_field_minima)
    : number_field_periods(input_field_periods),
      harmonic_part(input_field_periods, input_coefficients),
      order_zeros_field_maxima(input_order_zeros_field_maxima),
      order_zeros_field_minima(input_order_zeros_field_minima) {};

double general_curvature::zeros_func(double l) {
    double result = 1.0;
    for (int i = 0; i < order_zeros_field_maxima; i++) {
        result += std::sin(double(number_field_periods) * l / 2.0);
    }
    for (int i = 0; i < order_zeros_field_minima; i++) {
        result += std::cos(double(number_field_periods) * l / 2.0);
    }
    return result;
}

double general_curvature::half_period_function(double l) {
    return harmonic_part(l) * zeros_func(l);
};

double general_curvature::operator()(double l) {
    double result = 0.0;
    const double half_period_length =
        std::numbers::pi / double(number_field_periods);

    // map l into period
    while (l < 0) {
        l += 2 * half_period_length;
    }
    while (l > 2 * half_period_length) {
        l -= 2 * half_period_length;
    }

    // calculate function in first half period
    if (l <= half_period_length) {
        result += half_period_function(l);
    } else if (l > half_period_length) { // mirror function at half period
        result += half_period_function(half_period_length - l);
    }
    return result;
};

curve_system::curve_system(int input_field_periods,
                           std::vector<double> curvature_params,
                           std::vector<double> torsion_params,
                           int curvature_order_zeros_field_maxima,
                           int curvature_order_zeros_field_minima)
    : number_field_periods(input_field_periods),
      curvature(input_field_periods, curvature_params,
                curvature_order_zeros_field_maxima,
                curvature_order_zeros_field_minima),
      torsion(input_field_periods, torsion_params) {};
