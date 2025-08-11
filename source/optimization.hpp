#include "solver.hpp"
#include "system_class.hpp"
#include <vector>

double curve_closing_loss(std::vector<frenet_serret_frame> fs_solution);
double curve_closing_loss(curve_system sys);
