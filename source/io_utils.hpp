#include "system_class.hpp"
#include <string>
#include <vector>

void write_solution_to_file(std::string file_name,
                            std::vector<frenet_serret_frame> sample,
                            std::vector<double> arc_lengths);

void write_solution_to_file(std::string file_name, curve_system sys);
