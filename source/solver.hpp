#pragma once
constexpr int sample_size = 1000;
constexpr double low_integration_bound = 0.0;
constexpr double high_integration_bound = 2.0 * 3.141592653589;
constexpr double step_size = (high_integration_bound - low_integration_bound) /
                             static_cast<double>(sample_size);
constexpr int dims = 3;

#include <array>
#include <functional>
#include <vector>

typedef std::array<double, dims> real_vector;
typedef std::array<double, dims * dims> basis;
typedef std::array<double, (dims + 1) * dims> frenet_serret_frame;

// typedef for integrator
typedef frenet_serret_frame state_type;

void solve_frenet_serret(std::vector<state_type> &solution,
                         std::vector<double> &arc_lengths,
                         std::function<double(double)> curvature,
                         std::function<double(double)> torsion);

std::vector<real_vector> generate_curve(std::function<double(double)> curvature,
                                        std::function<double(double)> torsion);
