#pragma once
#include <functional>
#include <vector>

struct harmonic_series_pure_cos {
    int number_field_periods;
    std::vector<double> harmonic_coefficients;

    harmonic_series_pure_cos(int input_field_periods,
                             std::vector<double> input_coefficients);

    double operator()(double l);
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
};

typedef harmonic_series_pure_cos general_torsion;

struct curve_system {
    int number_field_periods;
    general_curvature curvature;
    general_torsion torsion;

    curve_system(int input_field_periods, std::vector<double> curvature_params,
                 std::vector<double> torsion_params,
                 int curvature_order_zeros_field_maxima,
                 int curvature_order_zeros_field_minima);
};
