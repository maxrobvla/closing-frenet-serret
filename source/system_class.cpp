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
                  std::cos(double(i) * double(number_field_periods) * l);
    }
    return result;
}

double
harmonic_series_pure_cos::derivative_wrt_coefficient(double l,
                                                     size_t coefficient_idx) {
    if (coefficient_idx == 0) {
        return 1.0;
    } else {
        return std::cos(double(coefficient_idx) * double(number_field_periods) *
                        l);
    }
}

bool check_harmonic_series_for_zeros(harmonic_series_pure_cos series) {
    double half_period_length =
        std::numbers::pi / double(series.number_field_periods);
    double DX = half_period_length / double(sample_size);
    double l;

    std::array<double, sample_size> grid_values;
    for (size_t i = 0; i < sample_size; i++) {
        l = DX * double(i);
        grid_values[i] = series(l);
    }

    // check the sign of every component
    // idea: if the sign changes throughout the half field period, there are
    // zeros inside the half field period
    bool has_zeros = true;
    bool positive = (grid_values[0] > 0);
    for (size_t i = 1; i < sample_size; i++) {
        if ((positive and (grid_values[i] < 0)) or
            ((not positive) and (grid_values[i] > 0))) {
            has_zeros = false;
            break;
        }
    }
    return has_zeros;
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
        result *= std::sin(double(number_field_periods) * l / 2.0);
    }
    for (int i = 0; i < order_zeros_field_minima; i++) {
        result *= std::cos(double(number_field_periods) * l / 2.0);
    }
    return result;
}

double general_curvature::half_period_function(double l) {
    return harmonic_part(l) * zeros_func(l);
}
double
general_curvature::derivative_half_period_function(double l,
                                                   size_t coefficient_idx) {
    return zeros_func(l) *
           harmonic_part.derivative_wrt_coefficient(l, coefficient_idx);
}

double general_curvature::derivative_wrt_coefficient(double l,
                                                     size_t coefficient_idx) {
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
        result += derivative_half_period_function(l, coefficient_idx);
    } else if (l > half_period_length) { // mirror function at half period
                                         // sign difference because of symmetry
        result += -derivative_half_period_function(half_period_length - l,
                                                   coefficient_idx);
    }
    return result;
}

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
      torsion(input_field_periods, torsion_params) {
    total_number_parameters = curvature_params.size() + torsion_params.size();
};

frenet_serret_frame curve_system::interpolate(double l) {
    {
        check_curve();
        frenet_serret_frame result;
        for (int i = 0; i < dims * (dims + 1); i++) {
            result[i] = interpolators[i](l);
        }
        return result;
    }
}

// set the actual curve and prepare interpolators
void curve_system::set_curve() {
    std::vector<frenet_serret_frame> fs_solution;
    std::vector<double> arc_lengths;
    solve_frenet_serret(fs_solution, arc_lengths, curvature, torsion);
    curve = std::move(fs_solution);

    std::vector<frenet_serret_frame> fs_derivatives(arc_lengths.size());
    for (size_t i = 0; i < arc_lengths.size(); i++) {
        frenet_serret_rhs(curve[i], fs_derivatives[i], arc_lengths[i],
                          curvature, torsion);
    }

    for (int d = 0; d < dims * (dims + 1); d++) {
        std::vector<std::array<double, 2>> coordinate_deriv_pair(
            arc_lengths.size());
        for (size_t i = 0; i < arc_lengths.size(); i++) {
            coordinate_deriv_pair[i] = {curve[i][d], fs_derivatives[i][d]};
        }
        interpolators.emplace_back(std::move(coordinate_deriv_pair),
                                   low_integration_bound, step_size);
    }
}

void curve_system::check_curve(bool verbose) {
    if ((curve.size() != sample_size + 1) or
        (interpolators.size() != dims * (dims + 1))) {
        if (verbose) {
            std::cout << "The curve system has not been set up correctly yet."
                      << std::endl;
        }
        set_curve();
    } else if (verbose) {
        std::cout << "The curve system has been set up correctly." << std::endl;
    }
}
