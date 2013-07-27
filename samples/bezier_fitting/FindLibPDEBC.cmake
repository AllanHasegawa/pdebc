
find_package(PkgConfig)
pkg_check_modules(PC_LIBPDEBC QUIET pdebc)
set(LIBPDEBC_DEFINITIONS ${PC_LIBPDEBC_CFLAGS_OTHER})

find_path(LIBPDEBC_INCLUDE_DIR pdebc/SequentialDE.hpp
          HINTS ${PC_LIBPDEBC_INCLUDEDIR} ${PC_LIBPDEBC_INCLUDE_DIRS}
          )

find_library(LIBPDEBC_LIBRARY NAMES pdebc
             HINTS ${PC_LIBPDEBC_LIBDIR} ${PC_LIBPDEBC_LIBRARY_DIRS}
             PATH_SUFFIXES pdebc )

set(LIBPDEBC_LIBRARIES "${LIBPDEBC_LIBRARY}")
set(LIBPDEBC_INCLUDE_DIRS ${LIBPDEBC_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LibPDEBC  DEFAULT_MSG
                                  LIBPDEBC_LIBRARY LIBPDEBC_INCLUDE_DIR)

mark_as_advanced(LIBPDEBC_INCLUDE_DIR LIBPDEBC_LIBRARY )