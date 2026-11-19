find_path(_boost_include_dir boost/version.hpp
  HINTS "${Boost_ROOT}" "$ENV{BOOST_ROOT}"
  PATH_SUFFIXES include)
if(NOT _boost_include_dir)
  set(PACKAGE_VERSION_COMPATIBLE FALSE)
  return()
endif()

file(STRINGS "${_boost_include_dir}/boost/version.hpp" _version_line
  REGEX "^#define BOOST_VERSION [0-9]+")
string(REGEX MATCH "[0-9]+" _version "${_version_line}")
math(EXPR _major "${_version} / 100000")
math(EXPR _minor "${_version} / 100 % 1000")
math(EXPR _patch "${_version} % 100")
set(PACKAGE_VERSION "${_major}.${_minor}.${_patch}")

if(PACKAGE_FIND_VERSION VERSION_GREATER PACKAGE_VERSION)
  set(PACKAGE_VERSION_COMPATIBLE FALSE)
else()
  set(PACKAGE_VERSION_COMPATIBLE TRUE)
  if(PACKAGE_FIND_VERSION VERSION_EQUAL PACKAGE_VERSION)
    set(PACKAGE_VERSION_EXACT TRUE)
  endif()
endif()
