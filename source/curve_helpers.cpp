#include "solver.hpp"
#include <cmath>

scalar_sample sample_function(double (*func)(double)) {
    scalar_sample result;
    for (int i = 0; i < sample_size; i++) {
        result[i] = func(step_size * i);
    }
    return result;
}

basis generate_starting_basis(double (*radius)(double),
                              double (*height)(double)) {
    // the curve is assumed to be parameterized in cylindrical coordinates,
    // where the angle is the free parameter
    // calculation is based on finite differences (to first order)
    // https://en.wikipedia.org/wiki/Finite_difference_coefficient
    basis starting_basis;

    // calculate tangent (first derivative of curve)
    real_vector tangent, normal, y0, y1, y2;
    y0 = {radius(0), 0, height(0)};
    y1 = {radius(step_size) * std::cos(step_size),
          radius(step_size) * std::sin(step_size), height(step_size)};
    y2 = {radius(2 * step_size) * std::cos(2 * step_size),
          radius(2 * step_size) * std::sin(2 * step_size),
          height(2 * step_size)};

    // don't need to devide by step_size because it will get normalized anyways
    for (int i = 0; i < 3; i++) {
        tangent[i] = y1[i] - y0[i];
    }
    // normalize tangent
    double tangent_length =
        std::sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1] +
                  tangent[2] * tangent[2]);
    for (int i = 0; i < 3; i++) {
        starting_basis[0][i] = tangent[i] / tangent_length;
    }

    // calculate normal (second derivative of curve)
    for (int i = 0; i < 3; i++) {
        normal[i] = y0[i] - 2 * y1[i] + y2[i];
    }
    double normal_length = std::sqrt(
        normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    for (int i = 0; i < 3; i++) {
        starting_basis[1][i] = normal[i] / normal_length;
    }

    // calculate binormal (cross product of tangent and normal)
    // cross product of two unit vectors is already normalized
    starting_basis[2][0] = starting_basis[0][1] * starting_basis[1][2] -
                           starting_basis[0][2] * starting_basis[1][1];
    starting_basis[2][1] = starting_basis[0][2] * starting_basis[1][0] -
                           starting_basis[0][0] * starting_basis[1][2];
    starting_basis[2][2] = starting_basis[0][0] * starting_basis[1][1] -
                           starting_basis[0][1] * starting_basis[1][0];

    return starting_basis;
}
