cmake_minimum_required(VERSION 3.20)

set(EACXX_VER 0.0.1)

project(
  eacxx
  VERSION ${EACXX_VER}
  LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
message(STATUS "Configuring {eacxx} version ${EACXX_VER}")

set(ROOT_DIR ${CMAKE_SOURCE_DIR})
set(SOURCE_DIR ${ROOT_DIR}/include)

# CMake configuration
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(compile_configuration)
include(GNUInstallDirs)

file(
  GLOB_RECURSE
  SOURCES
  CONFIGURE_DEPENDS
  ${SOURCE_DIR}/*.h
  ${SOURCE_DIR}/**/*.h
  ${SOURCE_DIR}/*.hpp
  ${SOURCE_DIR}/**/*.hpp)

add_library(eacxx INTERFACE)
target_include_directories(
  eacxx INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                  $<INSTALL_INTERFACE:include>)
set_target_properties(eacxx PROPERTIES LINKER_LANGUAGE CXX CXX_STANDARD 20)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/eacxxConfigVersion.cmake"
  VERSION ${EACXX_VER}
  COMPATIBILITY AnyNewerVersion)

install(
  TARGETS eacxx
  EXPORT eacxxTargets
  LIBRARY DESTINATION lib COMPONENT Runtime
  ARCHIVE DESTINATION lib COMPONENT Development
  RUNTIME DESTINATION bin COMPONENT Runtime
  PUBLIC_HEADER DESTINATION include COMPONENT Development
  BUNDLE DESTINATION bin COMPONENT Runtime)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/eacxxConfig.cmake.in"
  "${PROJECT_BINARY_DIR}/eacxxConfig.cmake" INSTALL_DESTINATION lib/cmake/mylib)

install(
  EXPORT eacxxTargets
  NAMESPACE eacxx::
  DESTINATION lib/cmake/eacxx)
install(FILES "${PROJECT_BINARY_DIR}/eacxxConfigVersion.cmake"
              "${PROJECT_BINARY_DIR}/eacxxConfig.cmake"
        DESTINATION lib/cmake/eacxx)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)
