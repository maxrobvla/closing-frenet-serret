#include "solver.hpp"
#include "system_class.hpp"
#include <cmath>
#include <vector>

double curve_closing_loss(std::vector<frenet_serret_frame> fs_solution) {
    double loss = 0.0;

    // get euclidian metric of position, tangent, normal
    // but allow for sign difference in normal

    // entries of the arrays
    // 0:position, 1:tangent, 2:normal, 3:normal inverted, 4: binormal, 5:
    // binormal inverted
    std::array<real_vector, 6> diff_vectors = {0.0};
    for (int i = 0; i < dims; i++) {
        diff_vectors[0][i] = fs_solution.front()[i] - fs_solution.back()[i];
        diff_vectors[1][i] =
            fs_solution.front()[dims + i] - fs_solution.back()[dims + i];
        diff_vectors[2][i] = fs_solution.front()[2 * dims + i] -
                             fs_solution.back()[2 * dims + i];
        diff_vectors[3][i] = fs_solution.front()[2 * dims + i] +
                             fs_solution.back()[2 * dims + i];
        diff_vectors[4][i] = fs_solution.front()[3 * dims + i] -
                             fs_solution.back()[3 * dims + i];
        diff_vectors[5][i] = fs_solution.front()[3 * dims + i] +
                             fs_solution.back()[3 * dims + i];
    }
    std::array<double, 6> single_losses = {0.0};
    for (int i = 0; i < dims; i++) {
        for (int j = 0; j < 6; j++) {
            single_losses[j] += diff_vectors[j][i] * diff_vectors[j][i];
        }
    }
    loss += std::sqrt(single_losses[0]);
    loss += 1e-1 * std::sqrt(single_losses[1]);
    if (single_losses[2] < single_losses[3]) {
        loss += 1e-2 * std::sqrt(single_losses[2]);
    } else {
        loss += 1e-2 * std ::sqrt(single_losses[3]);
    }
    if (single_losses[4] < single_losses[5]) {
        loss += 1e-2 * std::sqrt(single_losses[4]);
    } else {
        loss += 1e-2 * std::sqrt(single_losses[5]);
    }

    return loss;
}

// wrapper for curve_closing_loss
double curve_closing_loss(curve_system sys) {
    std::vector<frenet_serret_frame> solution;
    std::vector<double> arc_lengths;
    solve_frenet_serret(solution, arc_lengths, sys.curvature, sys.torsion);
    return curve_closing_loss(solution);
}
