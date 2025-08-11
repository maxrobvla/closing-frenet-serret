#include "argparse.hpp"
#include "io_utils.hpp"
#include "optimization.hpp"
#include "solver.hpp"
#include "system_class.hpp"
#include <cmath>
#include <iostream>

int main(int argc, char *argv[]) {

    // auto curvature = [](double) { return 1.0; };
    // auto torsion = [](double) { return 0.0; };
    // solve_frenet_serret(solution, arc_lengths, curvature, torsion);

    curve_system sys_eqs = parse_args(argc, argv);

    std::vector<frenet_serret_frame> solution;
    std::vector<double> arc_lengths;
    solve_frenet_serret(solution, arc_lengths, sys_eqs.curvature,
                        sys_eqs.torsion);

    write_solution_to_file("solution", solution, arc_lengths);

    std::cout << "Loss = " << curve_closing_loss(solution) << std::endl;
    return 0;
}
