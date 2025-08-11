#include "solver.hpp"
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/dense_output_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/generation/make_dense_output.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <functional>
#include <numbers>
#include <vector>

namespace ode = boost::numeric::odeint;

const double PI = std::numbers::pi;

void frenet_serret_rhs(const state_type &current_basis, state_type &derivative,
                       double l, std::function<double(double)> curvature,
                       std::function<double(double)> torsion) {
    for (int i = 0; i < 3; i++) {
        // position
        derivative[i] = current_basis[dims + i];
        // tangent
        derivative[dims + i] = curvature(l) * current_basis[2 * dims + i];
        // normal
        derivative[2 * dims + i] = -curvature(l) * current_basis[dims + i] +
                                   torsion(l) * current_basis[3 * dims + i];
        // binormal
        derivative[3 * dims + i] = -torsion(l) * current_basis[2 * dims + i];
    }
}

struct push_back_state_and_time {
    std::vector<state_type> &m_states;
    std::vector<double> &m_times;

    push_back_state_and_time(std::vector<state_type> &states,
                             std::vector<double> &times)
        : m_states(states), m_times(times) {}

    void operator()(const state_type &x, double t) {
        m_states.push_back(x);
        m_times.push_back(t);
    }
};

void solve_frenet_serret(std::vector<state_type> &solution,
                         std::vector<double> &arc_lengths,
                         std::function<double(double)> curvature,
                         std::function<double(double)> torsion) {
    state_type starting_frame;
    // 3dim. vectors: position, tangent, normal, binormal
    starting_frame = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};

    auto reduced_rhs = [curvature, torsion](const state_type &current_basis,
                                            state_type &derivative, double l) {
        frenet_serret_rhs(current_basis, derivative, l, curvature, torsion);
    };

    // setup integration method
    // ode::make_dense_output(absolute_tolerance, relative_tolerance,
    // integration_method)
    auto stepper = ode::make_dense_output(
        1e-9, 1e-9, ode::runge_kutta_dopri5<state_type>());

    size_t _ = ode::integrate_const(
        stepper, reduced_rhs, starting_frame, low_integration_bound,
        high_integration_bound, step_size,
        push_back_state_and_time(solution, arc_lengths));
}
