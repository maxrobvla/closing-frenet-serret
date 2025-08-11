#include "argparse.hpp"
#include "system_class.hpp"
#include <functional>
#include <iostream>
#include <string>
#include <vector>

curve_system parse_args(int argc, char *argv[]) {
    // parse inputs
    bool got_torsion = false;
    bool got_curvature = false;
    std::vector<double> torsion_params;
    std::vector<double> curvature_params;

    bool got_field_period = false;
    int field_period;

    bool got_order_zeros = false;
    int order_zeros_maxima;
    int order_zeros_minima;

    int active_param = 0;

    std::string arg;
    for (int i = 1; i < argc; i++) {
        arg = argv[i];
        if (arg == "-fp") {
            got_field_period = true;
            field_period = std::atoi(argv[i + 1]);
            i++; // skip next entry
        } else if (arg == "-c") {
            got_curvature = true;
            active_param = 1;
        } else if (arg == "-t") {
            got_torsion = true;
            active_param = 2;
        } else if (arg == "-zo") {
            got_order_zeros = true;
            order_zeros_maxima = std::atoi(argv[i + 1]);
            order_zeros_minima = std::atoi(argv[i + 2]);
            i += 2; // skip next two entries
        } else {
            switch (active_param) {
            case 1:
                curvature_params.push_back(std::atof(argv[i]));
                break;
            case 2:
                torsion_params.push_back(std::atof(argv[i]));
                break;
            default:
                std::cout << "Unexpected argument: " << argv[i] << std::endl;
            }
        }
    }

    // check if one entry is missing entirely
    if (not(got_field_period and got_curvature and got_torsion and
            got_order_zeros)) {
        std::exit(1.0);
    }

    // print read parameters to stdout for confirmation
    std::cout << "Field period " << field_period << std::endl;
    for (size_t i = 0; i < curvature_params.size(); i++) {
        std::cout << "Curvature " << i << "\t" << curvature_params[i]
                  << std::endl;
    }
    for (size_t i = 0; i < curvature_params.size(); i++) {
        std::cout << "Torsion " << i << "\t" << torsion_params[i] << std::endl;
    }

    // construct curve_sytem from inputs and return it
    curve_system sys_eqs(field_period, curvature_params, torsion_params,
                         order_zeros_maxima, order_zeros_minima);
    return sys_eqs;
}
