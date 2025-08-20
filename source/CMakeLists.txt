cmake_minimum_required(VERSION 3.22)
project(IPP-Frenet-Serret)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Enable warnings
if (MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Find OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    link_libraries(OpenMP::OpenMP_CXX)
endif()

set(CXX_SRC solver.cpp optimization.cpp io_utils.cpp system_class.cpp)
set(ARGP_SRC argparse.cpp)

# main executable
# add_executable("solver" main.cpp ${CXX_SRC} ${ARGP_SRC})

# build python module
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(closing_frenet_serret python_module.cpp ${CXX_SRC})
install(TARGETS closing_frenet_serret DESTINATION .)
