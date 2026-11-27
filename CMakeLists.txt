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
# Gate on the process exit code (cmd_validate returns 1 on any GS mismatch).
# Do not set PASS_REGULAR_EXPRESSION: CTest then ignores the exit code, so a
# future check that still prints "ground state correct for 17/17 elements" but
# returns non-zero (for example a lattice-constant failure) would be reported
# as a pass.
# "MISMATCH" is printed only on a wrong predicted ground state; it does not
# appear in a passing run's output.
set_tests_properties(validate PROPERTIES
  FAIL_REGULAR_EXPRESSION "MISMATCH")

if(ATOMFORGE_BUILD_TESTS)
  add_executable(test_atomforge tests/test_atomforge.cpp)
  target_link_libraries(test_atomforge PRIVATE atomforge)
  target_compile_options(test_atomforge PRIVATE -Wall -Wextra -O3)
  if(ATOMFORGE_NATIVE)
    target_compile_options(test_atomforge PRIVATE -march=native)
  endif()

  # One CTest entry per test group; labels allow e.g. `ctest -L derivatives` or `ctest -LE slow`.
  # The CLI test is named `validate` (from ci/actions-ctest) so the unit-test group is `validate_unit`.
  function(atomforge_add_group group)
    add_test(NAME ${group} COMMAND test_atomforge ${group})
    set_tests_properties(${group} PROPERTIES LABELS "${ARGN}" TIMEOUT 120)
  endfunction()
  atomforge_add_group(forces        unit derivatives)
  atomforge_add_group(virial        unit derivatives)
  atomforge_add_group(invariance    unit symmetry)
  atomforge_add_group(relax         physics)
  atomforge_add_group(born          physics)
  atomforge_add_group(hull          physics slow)
  atomforge_add_group(validate_unit unit)
  atomforge_add_group(cli            unit)
  atomforge_add_group(io             unit)

  # Regression tests for the NeighborList fix: atoms displaced by whole lattice vectors, and MD invariance under such shifts.
  atomforge_add_group(unwrapped     physics symmetry)
  atomforge_add_group(md_invariance physics symmetry)
endif()

# ExternalPotential bridge (include/atomforge/external.hpp), the native Lennard-Jones potential and quantum/worker.py.
# Groups that need python3 exit 77 when it is missing, which CTest reports as skipped.
if(ATOMFORGE_BUILD_TESTS)
  add_executable(test_external tests/test_external.cpp)
  target_link_libraries(test_external PRIVATE atomforge)
  target_compile_options(test_external PRIVATE -Wall -Wextra -O3)
  target_compile_definitions(test_external PRIVATE ATOMFORGE_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
  function(atomforge_add_external_group group)
    add_test(NAME ${group} COMMAND test_external ${group})
    set_tests_properties(${group} PROPERTIES LABELS "${ARGN}" TIMEOUT 300 SKIP_RETURN_CODE 77)
  endfunction()
  atomforge_add_external_group(lj              unit derivatives)
  atomforge_add_external_group(external        unit external)
  atomforge_add_external_group(external_relax  physics external)
  atomforge_add_external_group(external_errors unit external)
  add_test(NAME external_cli COMMAND test_external external_cli $<TARGET_FILE:atomforge_cli>)
  set_tests_properties(external_cli PROPERTIES LABELS "cli external" TIMEOUT 300 SKIP_RETURN_CODE 77)
endif()

# Finite-displacement phonons (include/atomforge/phonons.hpp).
if(ATOMFORGE_BUILD_TESTS)
  add_executable(test_phonons tests/test_phonons.cpp)
  target_link_libraries(test_phonons PRIVATE atomforge)
  target_compile_options(test_phonons PRIVATE -Wall -Wextra -O3)
  target_compile_definitions(test_phonons PRIVATE ATOMFORGE_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_test(NAME phonons COMMAND test_phonons phonons)
  set_tests_properties(phonons PROPERTIES LABELS "physics phonons" TIMEOUT 300)
  add_test(NAME phonons_lj COMMAND test_phonons phonons_lj)
  set_tests_properties(phonons_lj PROPERTIES LABELS "unit phonons" TIMEOUT 300)
endif()
