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

# set(SIM_SRC sim.cpp io_utils.cpp physics.cpp cell_lists.cpp particle.cpp parameters.cpp)  # General simulation source file
set(SIM_SRC solver.cpp)
set(ARGP_SRC argparse.cpp)

# Executable list: name, source file, compile definitions
set(EXE_LIST
  "solver" "${CMAKE_SOURCE_DIR}/main.cpp" "SOLVER=1"
)

set(OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/..)  # Output directory for executables

# Calculate number of executables
list(LENGTH EXE_LIST EXE_LIST_LEN)
math(EXPR EXE_COUNT "${EXE_LIST_LEN} / 3")

math(EXPR LOOP_RANGE "${EXE_COUNT}-1")

# Loop through the executables and set up each target
foreach(i RANGE 0 ${LOOP_RANGE})  # Fix: Use range from 0 to EXE_COUNT-1
    math(EXPR NAME_INDEX "${i} * 3")
    list(GET EXE_LIST ${NAME_INDEX} EXE_NAME)       # Get executable name
    math(EXPR SRC_INDEX "${i} * 3 + 1")
    list(GET EXE_LIST ${SRC_INDEX} SRC_FILE)     # Get source file
    math(EXPR DEF_INDEX "${i} * 3 + 2")
    list(GET EXE_LIST ${DEF_INDEX} TARGET_DEF)   # Get compile definitions

    # Add executable target
    add_executable(${EXE_NAME} ${SRC_FILE} ${SIM_SRC} ${ARGP_SRC})

    # Set macro definitions (convert space-separated to semicolon-separated list)
    string(REPLACE " " ";" TARGET_DEF_LIST ${TARGET_DEF})
    target_compile_definitions(${EXE_NAME} PUBLIC ${TARGET_DEF_LIST})

    # Add custom command to copy the executable after the build if it has changed
    # add_custom_command(
    #     TARGET ${EXE_NAME}
    #     POST_BUILD
    #     COMMAND ${CMAKE_COMMAND} -E copy_if_different
    #             $<TARGET_FILE:${EXE_NAME}> ${OUTPUT_DIRECTORY}
    #     COMMENT "Copying ${EXE_NAME} if it has changed"
    # )
endforeach()
