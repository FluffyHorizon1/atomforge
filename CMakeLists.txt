cmake_minimum_required(VERSION 3.16)
project(atomforge CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Host-specific codegen; leave OFF for reproducible CI and release binaries.
option(ATOMFORGE_NATIVE "Enable -march=native for the CLI" OFF)

add_library(atomforge INTERFACE)
target_include_directories(atomforge INTERFACE include)
add_executable(atomforge_cli src/main.cpp)
set_target_properties(atomforge_cli PROPERTIES OUTPUT_NAME atomforge)
target_link_libraries(atomforge_cli PRIVATE atomforge)
target_compile_options(atomforge_cli PRIVATE -Wall -Wextra -O3)
if(ATOMFORGE_NATIVE)
  target_compile_options(atomforge_cli PRIVATE -march=native)
endif()

enable_testing()
add_test(NAME validate COMMAND atomforge_cli validate)
# Gate on the process exit code (cmd_validate returns 1 on any GS mismatch).
# Do not set PASS_REGULAR_EXPRESSION: CTest then ignores the exit code, so a
# future check that still prints "ground state correct for 17/17 elements" but
# returns non-zero (for example a lattice-constant failure) would be reported
# as a pass.
# "MISMATCH" is printed only on a wrong predicted ground state; it does not
# appear in a passing run's output.
set_tests_properties(validate PROPERTIES
  FAIL_REGULAR_EXPRESSION "MISMATCH")
