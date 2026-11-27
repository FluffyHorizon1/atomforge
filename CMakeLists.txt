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

# --- forge.materials.v0 sidecar: record the git commit ("unknown" outside a git checkout).
find_package(Git QUIET)
set(ATOMFORGE_GIT_COMMIT "unknown")
if(GIT_FOUND)
  execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  OUTPUT_VARIABLE _af_commit OUTPUT_STRIP_TRAILING_WHITESPACE
                  RESULT_VARIABLE _af_rc ERROR_QUIET)
  if(_af_rc EQUAL 0 AND _af_commit)
    set(ATOMFORGE_GIT_COMMIT "${_af_commit}")
  endif()
endif()
# --- end forge.materials.v0 sidecar (commit) ---

add_library(atomforge INTERFACE)
target_include_directories(atomforge INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
# --- forge.materials.v0 sidecar: embed the commit in the sidecar writer.
target_compile_definitions(atomforge INTERFACE ATOMFORGE_GIT_COMMIT="${ATOMFORGE_GIT_COMMIT}")
# --- end forge.materials.v0 sidecar (commit define) ---
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
  atomforge_add_group(io             unit)

  # Regression tests for the NeighborList fix: atoms displaced by whole lattice vectors, and MD invariance under such shifts.
  atomforge_add_group(unwrapped     physics symmetry)
  atomforge_add_group(md_invariance physics symmetry)

  # Coverage groups: bookkeeping and error paths, file writers (POSCAR / extxyz / report.json), and the CLI itself.
  # `io` (registered above) now also runs the writer round trips.
  atomforge_add_group(structure_misc unit)
  # The cli group (a01 argument checks + the coverage checks) drives the real executable and runs a reduced `search`,
  # so it needs the path as a second argument.
  add_test(NAME cli COMMAND test_atomforge cli $<TARGET_FILE:atomforge_cli>)
  set_tests_properties(cli PROPERTIES LABELS "unit cli slow" TIMEOUT 600)
endif()

# --- forge.materials.v0 sidecar tests (self-contained block).
# Sits after the CLI validate test so it composes with ci/actions-ctest and with
# chore/atomforge-tests-suite: if both if(ATOMFORGE_BUILD_TESTS) blocks appear,
# keep both. Labels follow the suite's scheme (`ctest -L materials`, `-L unit`,
# `-LE slow`).
if(ATOMFORGE_BUILD_TESTS)
  add_executable(test_materials_contract tests/test_materials_contract.cpp)
  target_link_libraries(test_materials_contract PRIVATE atomforge)
  target_compile_options(test_materials_contract PRIVATE -Wall -Wextra -O3)
  if(ATOMFORGE_NATIVE)
    target_compile_options(test_materials_contract PRIVATE -march=native)
  endif()
  foreach(group schema placeholders determinism roundtrip)
    add_test(NAME materials_${group} COMMAND test_materials_contract ${group})
    set_tests_properties(materials_${group} PROPERTIES LABELS "materials;unit" TIMEOUT 120)
  endforeach()
  # Tiny Cu-Au search (no MC, 1 seed, short MD), run twice for determinism.
  add_test(NAME materials_real_search COMMAND test_materials_contract real_search)
  set_tests_properties(materials_real_search PROPERTIES LABELS "materials;physics;slow" TIMEOUT 300)
endif()
# --- end forge.materials.v0 sidecar tests ---
