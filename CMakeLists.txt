cmake_minimum_required(VERSION 3.22)
project(manifold_hull_repro LANGUAGES CXX)

include(CTest)

if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/manifold/CMakeLists.txt" OR
   NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/cgal/CGALConfig.cmake")
  message(FATAL_ERROR
    "Dependencies are missing. Run: git submodule update --init --recursive")
endif()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(MANIFOLD_TEST OFF CACHE BOOL "" FORCE)
set(MANIFOLD_CBIND OFF CACHE BOOL "" FORCE)
set(MANIFOLD_PYBIND OFF CACHE BOOL "" FORCE)
set(MANIFOLD_PAR OFF CACHE BOOL "" FORCE)
set(MANIFOLD_DOWNLOADS OFF CACHE BOOL "" FORCE)
set(ASSIMP_ENABLE OFF CACHE BOOL "" FORCE)
add_subdirectory(third_party/manifold manifold-build EXCLUDE_FROM_ALL)

set(CGAL_DISABLE_GMP ON CACHE BOOL "" FORCE)
set(Boost_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake" CACHE PATH "" FORCE)
find_package(CGAL 6.2 REQUIRED CONFIG
  PATHS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/cgal"
  NO_DEFAULT_PATH)

add_executable(manifold_hull_repro src/main.cpp)
target_link_libraries(manifold_hull_repro PRIVATE manifold CGAL::CGAL)

add_test(NAME manifold_hull_repro COMMAND manifold_hull_repro)
if(WIN32 AND MINGW AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.22)
  get_filename_component(_mingw_runtime_dir "${CMAKE_CXX_COMPILER}" DIRECTORY)
  set_tests_properties(manifold_hull_repro PROPERTIES
    ENVIRONMENT_MODIFICATION
      "PATH=path_list_prepend:${_mingw_runtime_dir}")
endif()


