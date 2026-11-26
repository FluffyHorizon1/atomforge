cmake_minimum_required(VERSION 3.16)
project(atomforge VERSION 0.1.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Host-specific codegen; leave OFF for reproducible CI and release binaries.
option(ATOMFORGE_NATIVE "Enable -march=native for the CLI" OFF)
option(ATOMFORGE_BUILD_TESTS "Build the tests/ CTest suite (the CLI validate test is always registered)" ON)

add_library(atomforge INTERFACE)
target_include_directories(atomforge INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_executable(atomforge_cli src/main.cpp)
set_target_properties(atomforge_cli PROPERTIES OUTPUT_NAME atomforge)
target_link_libraries(atomforge_cli PRIVATE atomforge)
target_compile_options(atomforge_cli PRIVATE -Wall -Wextra -O3)
if(ATOMFORGE_NATIVE)
  target_compile_options(atomforge_cli PRIVATE -march=native)
endif()

enable_testing()
add_test(NAME validate COMMAND atomforge_cli validate)
set_tests_properties(validate PROPERTIES
  PASS_REGULAR_EXPRESSION "ground state correct for 17/17 elements")
