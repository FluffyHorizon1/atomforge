cmake_minimum_required(VERSION 3.16)
project(atomforge CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
option(ATOMFORGE_BUILD_TESTS "Build the CTest suite under tests/" ON)
# Host-specific codegen; leave OFF for reproducible CI and release binaries.
option(ATOMFORGE_NATIVE "Enable -march=native for the CLI and tests" OFF)

# Record the git commit in forge.materials.v0 sidecars ("unknown" outside a git checkout).
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

add_library(atomforge INTERFACE)
target_include_directories(atomforge INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(atomforge INTERFACE ATOMFORGE_GIT_COMMIT="${ATOMFORGE_GIT_COMMIT}")
add_executable(atomforge_cli src/main.cpp)
set_target_properties(atomforge_cli PROPERTIES OUTPUT_NAME atomforge)
target_link_libraries(atomforge_cli PRIVATE atomforge)
target_compile_options(atomforge_cli PRIVATE -Wall -Wextra -O3)
if(ATOMFORGE_NATIVE)
  target_compile_options(atomforge_cli PRIVATE -march=native)
endif()

enable_testing()
add_test(NAME validate COMMAND atomforge_cli validate)
# Gate on the exit code (cmd_validate returns 1 on any GS mismatch); no PASS_REGULAR_EXPRESSION (it would override the exit code).
set_tests_properties(validate PROPERTIES FAIL_REGULAR_EXPRESSION "MISMATCH")

if(ATOMFORGE_BUILD_TESTS)
  add_executable(test_materials_contract tests/test_materials_contract.cpp)
  target_link_libraries(test_materials_contract PRIVATE atomforge)
  target_compile_options(test_materials_contract PRIVATE -Wall -Wextra -O3)
  if(ATOMFORGE_NATIVE)
    target_compile_options(test_materials_contract PRIVATE -march=native)
  endif()
  foreach(group schema placeholders determinism roundtrip real_search)
    add_test(NAME materials_${group} COMMAND test_materials_contract ${group})
    set_tests_properties(materials_${group} PROPERTIES LABELS "materials" TIMEOUT 300)
  endforeach()
endif()
