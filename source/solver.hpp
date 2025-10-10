#pragma once
constexpr unsigned int sample_size = 1000;
constexpr double low_integration_bound = 0.0f;
constexpr double high_integration_bound =
    2.0f * 3.141592653589793238462643383279502884L;
constexpr double step_size = (high_integration_bound - low_integration_bound) /
                             static_cast<double>(sample_size);
#include "system_class.hpp"
#include <functional>
#include <vector>

// typedef for integrator
typedef frenet_serret_frame state_type;

void frenet_serret_rhs(const state_type &current_basis, state_type &derivative,
                       double l, std::function<double(double)> curvature,
                       std::function<double(double)> torsion);

void solve_frenet_serret(std::vector<state_type> &solution,
                         std::vector<double> &arc_lengths,
                         std::function<double(double)> curvature,
                         std::function<double(double)> torsion);

void solve_frenet_serret(std::vector<state_type> &solution,
                         std::vector<double> &arc_lengths, curve_system sys);

#ifdef ENABLE_GRADIENT_CALCULATION
frenet_serret_frame calc_derivative_wrt_parameter(curve_system sys,
                                                  size_t parameter_idx);
#endif
