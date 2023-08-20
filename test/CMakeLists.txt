cmake_minimum_required(VERSION 3.20)

project(eacxx-tests LANGUAGES CXX)

find_package(eacxx REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(fmt 10.0.0 REQUIRED)
find_package(range-v3 REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main eacxx::eacxx fmt::fmt Eigen3::Eigen)
set_target_properties(main PROPERTIES LINKER_LANGUAGE CXX CXX_STANDARD 20)
target_compile_options(main PRIVATE -Wall -Wextra -O3 -march=native -ffast-math)
