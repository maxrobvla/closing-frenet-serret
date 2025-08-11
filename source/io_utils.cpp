#include "io_utils.hpp"
#include "solver.hpp"
#include "system_class.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void write_solution_to_file(std::string name,
                            std::vector<frenet_serret_frame> sample,
                            std::vector<double> arc_lengths) {
    std::filesystem::path save_path("data");
    std::string file_name = name + ".dat";
    std::ofstream out(save_path / file_name);
    out << "l,p1,p2,p3,t1,t2,t3,n1,n2,n3,b1,b2,b3" << std::endl;
    for (int i = 0; i < (int)arc_lengths.size(); i++) {
        out << arc_lengths[i] << ",";
        for (int j = 0; j < (3 * dims) + 2; j++) {
            out << sample[i][j] << ",";
        }
        out << sample[i][3 * dims + 2] << std::endl;
    }
}

void write_solution_to_file(std::string file_name, curve_system sys) {
    std::vector<frenet_serret_frame> solution;
    std::vector<double> arc_lengths;
    solve_frenet_serret(solution, arc_lengths, sys.curvature, sys.torsion);

    write_solution_to_file(file_name, solution, arc_lengths);
}
