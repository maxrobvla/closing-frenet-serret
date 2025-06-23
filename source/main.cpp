#include "solver.hpp"
#include <cmath>
#include <iostream>

int main(int argc, char *argv[]) {
    scalar_sample sampled_func = sample_function(std::cos);

    for (int i = 0; i < sample_size / 10; i++) {
        std::cout << sampled_func[i * 10] << std::endl;
    }
    return 0;
}
