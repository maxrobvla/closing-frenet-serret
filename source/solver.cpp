#include "solver.hpp"

basis fenet_serret_rhs(double curvature, double torsion, basis old_basis) {
    // https://en.wikipedia.org/wiki/Frenet-Serret_formulas
    basis result;
    for (int i = 0; i < 3; i++) {
        result[0][i] = curvature * old_basis[1][i];
        result[1][i] = -curvature * old_basis[0][i] + torsion * old_basis[2][i];
        result[2][i] = -torsion * old_basis[1][i];
    }
    return result;
}

basis_sample solve_fenet_serret(scalar_2x_sample curvature,
                                scalar_2x_sample torsion,
                                basis starting_basis) {
    basis_sample solution; // [time][vector][dimension]
    solution[0] = starting_basis;
    // TODO: finish the solver
    // https://en.wikipedia.org/wiki/Runge-Kutta_methods#The_Runge-Kutta_method
    basis k1, k2, k3, k4, yk;
    for (int s = 1; s < sample_size; s++) {
        // calculate all the helper terms
        k1 = fenet_serret_rhs(curvature[2 * (s - 1)], torsion[2 * (s - 1)],
                              solution[s]);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                yk[i][j] = solution[s - 1][i][j] + step_size / 2.0 * k1[i][j];
            }
        }
        k2 = fenet_serret_rhs(curvature[2 * s - 1], torsion[2 * s - 1], yk);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                yk[i][j] = solution[s - 1][i][j] + step_size / 2.0 * k2[i][j];
            }
        }
        k3 = fenet_serret_rhs(curvature[2 * s - 1], torsion[2 * s - 1], yk);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                yk[i][j] = solution[s - 1][i][j] + step_size * k3[i][j];
            }
        }
        k4 = fenet_serret_rhs(curvature[2 * s], torsion[2 * s], yk);

        // actual update
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                solution[s][i][j] =
                    solution[s - 1][i][j] +
                    step_size / 6.0 *
                        (k1[i][j] + 2.0 * (k2[i][j] + k3[i][j]) + k4[i][j]);
            }
        }
    }

    return solution;
}
