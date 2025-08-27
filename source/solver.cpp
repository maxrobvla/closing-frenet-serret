#include "solver.hpp"
#include "system_class.hpp"
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/symplectic_euler.hpp>
#include <functional>
#include <numbers>
#include <vector>

namespace ode = boost::numeric::odeint;

const double PI = std::numbers::pi;

void frenet_serret_rhs(const state_type &current_frame, state_type &derivative,
                       double l, std::function<double(double)> curvature,
                       std::function<double(double)> torsion) {
    for (int i = 0; i < 3; i++) {
        // position
        derivative[i] = current_frame[dims + i];
        // tangent
        derivative[dims + i] = curvature(l) * current_frame[2 * dims + i];
        // normal
        derivative[2 * dims + i] = -curvature(l) * current_frame[dims + i] +
                                   torsion(l) * current_frame[3 * dims + i];
        // binormal
        derivative[3 * dims + i] = -torsion(l) * current_frame[2 * dims + i];
    }
}

void frenet_serret_rhs(const state_type &current_frame, state_type &derivative,
                       double l, curve_system sys) {
    frenet_serret_rhs(current_frame, derivative, l, sys.curvature, sys.torsion);
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
                         std::vector<double> &arc_lengths, curve_system sys) {
    solve_frenet_serret(solution, arc_lengths, sys.curvature, sys.torsion);
}

void solve_frenet_serret(std::vector<state_type> &solution,
                         std::vector<double> &arc_lengths,
                         std::function<double(double)> curvature,
                         std::function<double(double)> torsion) {
    state_type starting_frame;
    // 3dim. vectors: position, tangent, normal, binormal
    starting_frame = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};

    auto reduced_rhs = [curvature, torsion](const state_type &current_frame,
                                            state_type &derivative, double l) {
        frenet_serret_rhs(current_frame, derivative, l, curvature, torsion);
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

void jacobian_frenet_serret_rhs(const state_type &frame_derivative,
                                state_type &derivative, double l,
                                curve_system sys, size_t parameter_idx) {
    // product rule derivative of frame
    frenet_serret_rhs(frame_derivative, derivative, l, sys);

    // product rule derivative of matrix
    if (parameter_idx > sys.total_number_parameters) {
        std::cerr << "Invalid parameter index when calculating derivative"
                  << std::endl;
        return;
    }
    frenet_serret_frame curve_at_point = sys.interpolate(l);
    if (parameter_idx <
        sys.curvature.harmonic_part.harmonic_coefficients.size()) {
        // add derivative with respect to curvature parameter
        for (int i = 0; i < 3; i++) {
            int coefficient_idx = parameter_idx;
            // tangent
            derivative[dims + i] +=
                sys.curvature.derivative_wrt_coefficient(l, coefficient_idx) *
                curve_at_point[2 * dims + i];
            // normal
            derivative[2 * dims + i] +=
                -sys.curvature.derivative_wrt_coefficient(l, coefficient_idx) *
                curve_at_point[dims + i];
        }
    } else {
        // add derivative with respect to torsion parameter
        int coefficient_idx =
            parameter_idx -
            sys.curvature.harmonic_part.harmonic_coefficients.size();
        for (int i = 0; i < 3; i++) {
            // normal
            derivative[2 * dims + i] +=
                sys.torsion.derivative_wrt_coefficient(l, coefficient_idx) *
                curve_at_point[3 * dims + i];
            // binormal
            derivative[3 * dims + i] +=
                -sys.torsion.derivative_wrt_coefficient(l, coefficient_idx) *
                curve_at_point[2 * dims + i];
        }
    }
}

frenet_serret_frame calc_derivative_wrt_parameter(curve_system sys,
                                                  size_t parameter_idx) {

    sys.check_curve();

    state_type starting_frame;
    // 3dim. vectors: position, tangent, normal, binormal
    starting_frame = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    auto reduced_rhs = [sys, parameter_idx](const state_type &current_frame,
                                            state_type &derivative, double l) {
        jacobian_frenet_serret_rhs(current_frame, derivative, l, sys,
                                   parameter_idx);
    };

    // setup integration method
    // ode::make_dense_output(absolute_tolerance, relative_tolerance,
    // integration_method)
    auto stepper = ode::make_dense_output(
        1e-9, 1e-9, ode::runge_kutta_dopri5<state_type>());

    std::vector<state_type> state;
    std::vector<double> time;

    size_t _ = ode::integrate_adaptive(stepper, reduced_rhs, starting_frame,
                                       low_integration_bound,
                                       high_integration_bound, step_size,
                                       push_back_state_and_time(state, time));
    return state.back();
}
