#include "solver.hpp"
#include "system_class.hpp"
#include <vector>

double curve_closing_loss(std::vector<frenet_serret_frame> fs_solution);
double curve_closing_loss(curve_system sys);

std::vector<double> jacobian_curve_closing_loss(curve_system sys);

std::vector<double>
curve_closing_loss_vector(std::vector<frenet_serret_frame> fs_solution);
std::vector<double> curve_closing_loss_vector(curve_system sys);

std::vector<std::vector<double>>
jacobian_curve_closing_loss_vector(curve_system sys,
                                   std::vector<frenet_serret_frame> curve);
std::vector<std::vector<double>>
jacobian_curve_closing_loss_vector(curve_system sys);
