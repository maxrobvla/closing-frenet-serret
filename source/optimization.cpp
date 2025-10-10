#include "optimization.hpp"
#include "solver.hpp"
#include "system_class.hpp"
#include <vector>

std::vector<double>
curve_closing_loss_vector(std::vector<frenet_serret_frame> fs_solution) {

    // dims*position, dims*tangent, 1*normal, 1*binormal
    std::vector<double> loss(2 * dims + 2, 0.0);

    // get euclidian metric of position, tangent end - start
    // but allow for sign difference in normal
    for (unsigned int d = 0; d < dims; d++) {
        loss[POSITION + d] = fs_solution.back()[POSITION + d] -
                             fs_solution.front()[POSITION + d];
        loss[TANGENT + d] =
            fs_solution.back()[TANGENT + d] - fs_solution.front()[TANGENT + d];
    }

    // 1 - (dot product (bi)normal)^2
    double normal_dot = 0.0;
    double binormal_dot = 0.0;

    for (int d = 0; d < dims; d++) {
        normal_dot +=
            fs_solution.back()[NORMAL + d] * fs_solution.front()[NORMAL + d];
        binormal_dot += fs_solution.back()[BINORMAL + d] *
                        fs_solution.front()[BINORMAL + d];
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
        loss += loss_vector[POSITION + i] * loss_vector[POSITION + i];
    }
    for (int i = 0; i < dims; i++) {
        loss += loss_vector[TANGENT + i] * loss_vector[TANGENT + i];
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

double sys_dot_product(curve_system sys, unsigned int component_1,
                       unsigned int component_2) {
    double result = 0;
    for (size_t d = 0; d < dims; d++) {
        result += sys.curve.front()[component_1 + d] *
                  sys.curve.back()[component_2 + d];
    }
    return result;
}

// different loss functions for use in root find
std::vector<double> curve_closing_loss_root(curve_system sys,
                                            unsigned int number_conditions) {
    sys.check_curve();

    std::vector<double> result;

    if (number_conditions == 0) {
        if (sys.number_field_periods == 1) {
            number_conditions = 1;
        } else if (sys.number_field_periods >= 2) {
            number_conditions = 2;
        }
    }

    result.reserve(number_conditions);

    // check if curve is closed (-2...2)
    if (number_conditions >= 3) {
        double diff = 0.0;
        for (size_t d = 0; d < dims; d++) {
            diff = sys.curve.front()[d] - sys.curve.back()[d];
            if (d == 0) {
                result.push_back(diff * diff);
            } else {
                result.back() += diff * diff;
            }
        }
    }
    // (0...2)
    if (number_conditions >= 4) {
        result.push_back(1 - sys_dot_product(sys, TANGENT, TANGENT));
    }

    // check if frames at start and end points are aligned (-1...1)
    if (number_conditions >= 1) {
        result.push_back(sys_dot_product(sys, TANGENT, NORMAL));
    }
    if (number_conditions >= 2) {
        result.push_back(sys_dot_product(sys, TANGENT, BINORMAL));
    }
    if (number_conditions >= 5) {
        result.push_back(sys_dot_product(sys, NORMAL, BINORMAL));
    }
    return result;
}

#ifdef ENABLE_GRADIENT_CALCULATION
// function for calculating jacobian of loss
std::vector<double> jacobian_curve_closing_loss(curve_system sys) {
    sys.check_curve(true);

    std::vector<double> result(sys.total_number_parameters, 0.0);

    std::vector<std::vector<double>> jacobian_loss_vector =
        jacobian_curve_closing_loss_vector(sys);

    size_t entry = 0.0;
    for (size_t i = 0; i < sys.total_number_parameters; i++) {
        for (int d = 0; d < dims; d++) {
            // position
            entry = POSITION + d;
            result[i] += 2.f *
                         (sys.curve.back()[entry] - sys.curve.front()[entry]) *
                         jacobian_loss_vector[i][entry];

            // tangent
            entry = TANGENT + d;
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

std::vector<std::vector<double>>
jacobian_curve_closing_loss_vector(curve_system sys) {
    sys.check_curve(true);

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
            entry = NORMAL + d;
            normal_dot_product +=
                sys.curve.back()[entry] * starting_frame[entry];

            entry = BINORMAL + d;
            binormal_dot_product +=
                sys.curve.back()[entry] * starting_frame[entry];
        }

        for (int d = 0; d < dims; d++) {
            // position
            entry = POSITION + d;
            result[i][entry] += derivative[entry];

            // tangent
            entry = TANGENT + d;
            result[i][entry] += derivative[entry];

            // normal
            entry = NORMAL + d;
            result[i][2 * dims] += -2.f * normal_dot_product *
                                   derivative[entry] * starting_frame[entry];

            // binormal
            entry = BINORMAL + d;
            result[i][2 * dims + 1] += -2.f * binormal_dot_product *
                                       derivative[entry] *
                                       starting_frame[entry];
        }
    }
    return result;
}
#endif
