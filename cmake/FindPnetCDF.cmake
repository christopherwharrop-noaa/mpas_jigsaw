# FindPnetCDF
# -----------
# Find PnetCDF (Parallel-NetCDF) library.
#
# Result variables:
#   PnetCDF_FOUND        - True if PnetCDF is found
#   PnetCDF_INCLUDE_DIRS - Include directories
#   PnetCDF_LIBRARIES    - Libraries to link
#   PnetCDF_VERSION      - Version string
#
# Imported targets:
#   PnetCDF::PnetCDF     - Imported target for PnetCDF
#
# Hints:
#   PNETCDF              - Environment variable pointing to install prefix
#   CMAKE_PREFIX_PATH    - Standard CMake search path

find_path(PnetCDF_INCLUDE_DIR pnetcdf.h
  HINTS ENV PNETCDF
  PATH_SUFFIXES include
)

find_library(PnetCDF_LIBRARY pnetcdf
  HINTS ENV PNETCDF
  PATH_SUFFIXES lib lib64
)

# Try to get version from pnetcdf.h
if(PnetCDF_INCLUDE_DIR AND EXISTS "${PnetCDF_INCLUDE_DIR}/pnetcdf.h")
  file(STRINGS "${PnetCDF_INCLUDE_DIR}/pnetcdf.h" _pnetcdf_version_line
    REGEX "^#define[ \t]+PNETCDF_VERSION[ \t]+\"[^\"]*\"")
  if(_pnetcdf_version_line)
    string(REGEX REPLACE ".*\"([^\"]*)\".*" "\\1" PnetCDF_VERSION "${_pnetcdf_version_line}")
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PnetCDF
  REQUIRED_VARS PnetCDF_LIBRARY PnetCDF_INCLUDE_DIR
  VERSION_VAR PnetCDF_VERSION
)

if(PnetCDF_FOUND AND NOT TARGET PnetCDF::PnetCDF)
  add_library(PnetCDF::PnetCDF UNKNOWN IMPORTED)
  set_target_properties(PnetCDF::PnetCDF PROPERTIES
    IMPORTED_LOCATION "${PnetCDF_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${PnetCDF_INCLUDE_DIR}"
  )
endif()

mark_as_advanced(PnetCDF_INCLUDE_DIR PnetCDF_LIBRARY)
