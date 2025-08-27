#include "solver.hpp"
#include "system_class.hpp"
#include <vector>

std::vector<double>
curve_closing_loss_vector(std::vector<frenet_serret_frame> fs_solution) {

    // dims*position, dims*tangent, 1*normal, 1*binormal
    std::vector<double> loss(2 * dims + 2, 0.0);

    // get euclidian metric of position, tangent end - start
    // but allow for sign difference in normal
    for (int i = 0; i < dims; i++) {
        loss[i] = fs_solution.back()[i] - fs_solution.front()[i];
        loss[dims + i] =
            fs_solution.back()[dims + i] - fs_solution.front()[dims + i];
    }

    // 1 - (dot product (bi)normal)^2
    double normal_dot = 0.0;
    double binormal_dot = 0.0;

    for (int i = 0; i < dims; i++) {
        normal_dot += fs_solution.back()[2 * dims + i] *
                      fs_solution.front()[2 * dims + i];
        binormal_dot += fs_solution.back()[3 * dims + i] *
                        fs_solution.front()[3 * dims + i];
    }
    loss[2 * dims] = 1 - normal_dot * normal_dot;
    loss[2 * dims + 1] = 1 - binormal_dot * binormal_dot;

    return loss;
}

double curve_closing_loss(std::vector<frenet_serret_frame> fs_solution) {
    double loss = 0.0;
    std::vector<double> loss_vector = curve_closing_loss_vector(fs_solution);

    // get metric squared of position, tangent (end - start)
    for (int i = 0; i < dims; i++) {
        loss += loss_vector[i] * loss_vector[i];
    }
    for (int i = 0; i < dims; i++) {
        loss += loss_vector[dims + i] * loss_vector[dims + i];
    }
    // 1 - (dot product (bi)normal)^2
    loss += loss_vector[2 * dims];
    loss += loss_vector[2 * dims + 1];

    return loss;
}

// wrapper for curve_closing_loss
double curve_closing_loss(curve_system sys) {
    sys.check_curve();
    return curve_closing_loss(sys.curve);
}

std::vector<double> curve_closing_loss_vector(curve_system sys) {
    sys.check_curve();
    return curve_closing_loss_vector(sys.curve);
}

std::vector<std::vector<double>>
jacobian_curve_closing_loss_vector(curve_system sys) {
    sys.check_curve();

    std::vector<std::vector<double>> result(
        sys.total_number_parameters, std::vector<double>(2 * dims + 2, 0.0));

    frenet_serret_frame derivative;
    frenet_serret_frame starting_frame = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
    size_t entry = 0.0;
    double normal_dot_product;
    double binormal_dot_product;
    for (size_t i = 0; i < sys.total_number_parameters; i++) {
        derivative = calc_derivative_wrt_parameter(sys, i);

        normal_dot_product = 0.0;
        binormal_dot_product = 0.0;
        for (int d = 0; d < dims; d++) {
            entry = 2 * dims + d;
            normal_dot_product +=
                sys.curve.back()[entry] * starting_frame[entry];

            entry = 3 * dims + d;
            binormal_dot_product +=
                sys.curve.back()[entry] * starting_frame[entry];
        }

        for (int d = 0; d < dims; d++) {
            // position
            entry = d;
            result[i][entry] += derivative[entry];

            // tangent
            entry = dims + d;
            result[i][entry] += derivative[entry];

            // normal
            entry = 2 * dims + d;
            result[i][2 * dims] += -2.f * normal_dot_product *
                                   derivative[entry] * starting_frame[entry];

            // binormal
            entry = 3 * dims + d;
            result[i][2 * dims + 1] += -2.f * binormal_dot_product *
                                       derivative[entry] *
                                       starting_frame[entry];
        }
    }
    return result;
}

// function for calculating jacobian of loss
// TODO: test using finite differences
std::vector<double> jacobian_curve_closing_loss(curve_system sys) {
    sys.check_curve();

    std::vector<double> result(sys.total_number_parameters, 0.0);

    std::vector<std::vector<double>> jacobian_loss_vector =
        jacobian_curve_closing_loss_vector(sys);

    size_t entry = 0.0;
    for (size_t i = 0; i < sys.total_number_parameters; i++) {
        for (int d = 0; d < dims; d++) {
            // position
            entry = d;
            result[i] += 2.f *
                         (sys.curve.back()[entry] - sys.curve.front()[entry]) *
                         jacobian_loss_vector[i][entry];

            // tangent
            entry = dims + d;
            result[i] += 2.f *
                         (sys.curve.back()[entry] - sys.curve.front()[entry]) *
                         jacobian_loss_vector[i][entry];
        }
        // normal dot product
        result[i] += jacobian_loss_vector[i][2 * dims];

        // binormal dot product
        result[i] += jacobian_loss_vector[i][2 * dims + 1];
    }
    return result;
}
