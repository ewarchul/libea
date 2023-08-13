cmake_minimum_required(VERSION 3.20)

set(LIBEACPP_VER 0.0.1)

project(
  eacpp
  VERSION ${LIBEACPP_VER}
  LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
message(STATUS "Configuring {libeacpp} version ${LIBEACPP_VER}")

set(ROOT_DIR ${CMAKE_SOURCE_DIR})
set(SOURCE_DIR ${ROOT_DIR}/include)

# Artifatcs directories
set(ARTIFACT_DIRECTORY
    ${CMAKE_SOURCE_DIR}/artifacts
    CACHE STRING "artifacts")
set(CMAKE_BINARY_DIR ${ARTIFACT_DIRECTORY}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/exe)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/archive)

# CMake configuration
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(compile_configuration)

# {libeacpp} dependencies
find_package(spdlog REQUIRED)
find_package(fmt 10.0.0 REQUIRED)
find_package(range-v3 REQUIRED)
find_package(Eigen3 3.0.0 REQUIRED)
find_package(blaze REQUIRED)
find_package(libcmaes REQUIRED)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package (Eigen3 3.4 REQUIRED NO_MODULE)

# file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS ${SOURCE_DIR}/*.cpp
# ${SOURCE_DIR}/**/*.cpp)
file(
  GLOB_RECURSE
  SOURCES
  CONFIGURE_DEPENDS
  ${SOURCE_DIR}/*.h
  ${SOURCE_DIR}/**/*.h
  ${SOURCE_DIR}/*.hpp
  ${SOURCE_DIR}/**/*.hpp)

add_library(eacpp INTERFACE)
target_include_directories(eacpp INTERFACE include)
# add_library(eacpp SHARED ${SOURCES})
set_target_properties(eacpp PROPERTIES LINKER_LANGUAGE CXX CXX_STANDARD 20)

add_executable(main main.cpp)
target_link_libraries(main eacpp spdlog::spdlog fmt::fmt BLAS::BLAS LAPACK::LAPACK
  blaze::blaze Eigen3::Eigen)
set_target_properties(main PROPERTIES CXX_STANDARD 20)
target_compile_options(main PRIVATE ${PEDANTIC_COMPILE_FLAGS}
                                    ${OPTIM_COMPILE_FLAGS})
