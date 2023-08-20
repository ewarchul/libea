cmake_minimum_required(VERSION 3.16)

project(cmaes_benchmark LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(eacxx REQUIRED)
find_package(libcmaes REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(fmt 10.0.0 REQUIRED)
find_package(range-v3 REQUIRED)


add_executable(cmaes_benchmark benchmark.cpp)
target_link_libraries(cmaes_benchmark eacxx::eacxx libcmaes::cmaes Eigen3::Eigen fmt::fmt)
set_target_properties(cmaes_benchmark PROPERTIES CXX_STANDARD 20)
target_compile_options(cmaes_benchmark PRIVATE -O3 -ffast-math -march=native )


