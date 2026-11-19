# Minimal header-only Boost package for CGAL. CMake 4 no longer ships
# FindBoost, while this reproducer needs no compiled Boost libraries.
find_path(Boost_INCLUDE_DIR boost/version.hpp
  HINTS "${Boost_ROOT}" "$ENV{BOOST_ROOT}"
  PATH_SUFFIXES include)

if(NOT Boost_INCLUDE_DIR)
  set(Boost_FOUND FALSE)
  set(Boost_NOT_FOUND_MESSAGE
    "Boost headers were not found; install Boost or set BOOST_ROOT")
  return()
endif()

set(Boost_INCLUDE_DIRS "${Boost_INCLUDE_DIR}")
set(Boost_LIBRARIES "")
foreach(_target Boost::headers Boost::boost)
  if(NOT TARGET ${_target})
    add_library(${_target} INTERFACE IMPORTED)
    set_target_properties(${_target} PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIR}")
  endif()
endforeach()
unset(_target)
set(Boost_FOUND TRUE)
