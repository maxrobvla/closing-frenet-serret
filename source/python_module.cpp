#include "io_utils.hpp"
#include "optimization.hpp"
#include "solver.hpp"
#include "system_class.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

int curvature_harmonic_part_check_zeros(
    int number_field_periods, py::array_t<double> curvature_coefficients) {

    std::vector<double> curvature_coefficients_vector,
        torsion_coefficients_vector;

    for (int i = 0; i < curvature_coefficients.size(); i++) {
        curvature_coefficients_vector.push_back(curvature_coefficients.at(i));
    }

    return check_harmonic_series_for_zeros(number_field_periods,
                                           curvature_coefficients_vector);
}

curve_system
create_system_from_python(int number_field_periods,
                          py::array_t<int> order_zeros,
                          py::array_t<double> curvature_coefficients,
                          py::array_t<double> torsion_coefficients) {
    std::vector<double> curvature_coefficients_vector,
        torsion_coefficients_vector;

    for (int i = 0; i < curvature_coefficients.size(); i++) {
        curvature_coefficients_vector.push_back(curvature_coefficients.at(i));
    }

    for (int i = 0; i < torsion_coefficients.size(); i++) {
        torsion_coefficients_vector.push_back(torsion_coefficients.at(i));
    }

    curve_system sys(number_field_periods, curvature_coefficients_vector,
                     torsion_coefficients_vector, order_zeros.at(0),
                     order_zeros.at(1));
    return sys;
}

double python_loss(int number_field_periods, py::array_t<int> order_zeros,
                   py::array_t<double> curvature_coefficients,
                   py::array_t<double> torsion_coefficients) {
    curve_system sys =
        create_system_from_python(number_field_periods, order_zeros,
                                  curvature_coefficients, torsion_coefficients);

    return curve_closing_loss(sys);
}

void save_curve(std::string filename, int number_field_periods,
                py::array_t<int> order_zeros,
                py::array_t<double> curvature_coefficients,
                py::array_t<double> torsion_coefficients) {
    curve_system sys =
        create_system_from_python(number_field_periods, order_zeros,
                                  curvature_coefficients, torsion_coefficients);

    write_solution_to_file(filename, sys);
}

std::vector<double>
python_loss_vector(int number_field_periods, py::array_t<int> order_zeros,
                   py::array_t<double> curvature_coefficients,
                   py::array_t<double> torsion_coefficients) {
    curve_system sys =
        create_system_from_python(number_field_periods, order_zeros,
                                  curvature_coefficients, torsion_coefficients);

    return curve_closing_loss_vector(sys);
}

py::tuple python_get_curve(int number_field_periods,
                           py::array_t<int> order_zeros,
                           py::array_t<double> curvature_coefficients,
                           py::array_t<double> torsion_coefficients) {
    curve_system sys =
        create_system_from_python(number_field_periods, order_zeros,
                                  curvature_coefficients, torsion_coefficients);
    sys.check_curve();
    py::array_t<double> curve_x(sys.curve.size());
    py::array_t<double> curve_y(sys.curve.size());
    py::array_t<double> curve_z(sys.curve.size());
    for (size_t i = 0; i < sys.curve.size(); i++) {
        curve_x.mutable_at(i) = sys.curve[i][0];
        curve_y.mutable_at(i) = sys.curve[i][1];
        curve_z.mutable_at(i) = sys.curve[i][2];
    }
    return py::make_tuple(curve_x, curve_y, curve_z);
}

std::vector<double>
python_loss_root_find(int number_field_periods, py::array_t<int> order_zeros,
                      py::array_t<double> curvature_coefficients,
                      py::array_t<double> torsion_coefficients,
                      unsigned int number_conditions) {
    curve_system sys =
        create_system_from_python(number_field_periods, order_zeros,
                                  curvature_coefficients, torsion_coefficients);
    sys.check_curve();

    return curve_closing_loss_root(sys, number_conditions);
}

#ifdef ENABLE_GRADIENT_CALCULATION
py::array_t<double>
python_jacobian_of_loss(int number_field_periods, py::array_t<int> order_zeros,
                        py::array_t<double> curvature_coefficients,
                        py::array_t<double> torsion_coefficients) {
    curve_system sys =
        create_system_from_python(number_field_periods, order_zeros,
                                  curvature_coefficients, torsion_coefficients);

    std::vector<frenet_serret_frame> curve;
    std::vector<double> arc_lengths;

    solve_frenet_serret(curve, arc_lengths, sys);

    std::vector<double> jacobian_vec = jacobian_curve_closing_loss(sys);

    auto jacobian = py::array_t<double>(jacobian_vec.size());
    auto jacobian_buffer = jacobian.request();
    auto jacobian_ptr = (double *)jacobian_buffer.ptr;

    std::memcpy((double *)jacobian_ptr, jacobian_vec.data(),
                jacobian_vec.size() * sizeof(double));

    return jacobian;
}

std::vector<std::vector<double>>
python_jacobian_of_loss_vector(int number_field_periods,
                               py::array_t<int> order_zeros,
                               py::array_t<double> curvature_coefficients,
                               py::array_t<double> torsion_coefficients) {
    curve_system sys =
        create_system_from_python(number_field_periods, order_zeros,
                                  curvature_coefficients, torsion_coefficients);

    return jacobian_curve_closing_loss_vector(sys);
}

py::tuple
python_jacobian_and_scalar_loss(int number_field_periods,
                                py::array_t<int> order_zeros,
                                py::array_t<double> curvature_coefficients,
                                py::array_t<double> torsion_coefficients) {

    curve_system sys =
        create_system_from_python(number_field_periods, order_zeros,
                                  curvature_coefficients, torsion_coefficients);
    sys.check_curve();
    double loss = curve_closing_loss(sys);
    std::vector<double> jacobian = jacobian_curve_closing_loss(sys);
    return py::make_tuple(loss, jacobian);
}

py::tuple
python_jacobian_and_vector_loss(int number_field_periods,
                                py::array_t<int> order_zeros,
                                py::array_t<double> curvature_coefficients,
                                py::array_t<double> torsion_coefficients) {

    curve_system sys =
        create_system_from_python(number_field_periods, order_zeros,
                                  curvature_coefficients, torsion_coefficients);
    sys.check_curve();
    std::vector<double> loss = curve_closing_loss_vector(sys);
    std::vector<std::vector<double>> jacobian =
        jacobian_curve_closing_loss_vector(sys);
    return py::make_tuple(loss, jacobian);
}
#endif

PYBIND11_MODULE(closing_frenet_serret, m) {
    m.def("loss", python_loss,
          "Returns curve closing loss for curvature and torsion specified by "
          "given parameters",
          py::arg("number_field_periods"), py::arg("order_zeros"),
          py::arg("curvature_coefficients"), py::arg("torsion_coefficients"));
    m.def("save_curve", save_curve,
          "Saves curve generated by curvature and torsion specified by given "
          "parameters to a file with name 'data/{filename}.dat'",
          py::arg("filename"), py::arg("number_field_periods"),
          py::arg("order_zeros"), py::arg("curvature_coefficients"),
          py::arg("torsion_coefficients"));
    m.def("curvature_harmonic_part_check_zeros",
          curvature_harmonic_part_check_zeros, py::arg("number_field_periods"),
          py::arg("curvature_coefficients"));
    m.def("loss_vector", python_loss_vector,
          "Returns curve closing loss for curvature and torsion specified by "
          "given parameters",
          py::arg("number_field_periods"), py::arg("order_zeros"),
          py::arg("curvature_coefficients"), py::arg("torsion_coefficients"));
    m.def("get_curve", python_get_curve, py::arg("number_field_periods"),
          py::arg("order_zeros"), py::arg("curvature_coefficients"),
          py::arg("torsion_coefficients"));
    m.def("loss_root_find", python_loss_root_find,
          py::arg("number_field_periods"), py::arg("order_zeros"),
          py::arg("curvature_coefficients"), py::arg("torsion_coefficients"),
          py::arg("number_conditions"));
#ifdef ENABLE_GRADIENT_CALCULATION
    m.def("loss_jacobian", python_jacobian_of_loss,
          py::arg("number_field_periods"), py::arg("order_zeros"),
          py::arg("curvature_coefficients"), py::arg("torsion_coefficients"));
    m.def("loss_jacobian_vector", python_jacobian_of_loss_vector,
          py::arg("number_field_periods"), py::arg("order_zeros"),
          py::arg("curvature_coefficients"), py::arg("torsion_coefficients"));
    m.def("jacobian_and_scalar_loss", python_jacobian_and_scalar_loss,
          py::arg("number_field_periods"), py::arg("order_zeros"),
          py::arg("curvature_coefficients"), py::arg("torsion_coefficients"));
#endif
}
