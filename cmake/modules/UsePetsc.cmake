find_package(PkgConfig REQUIRED)
pkg_search_module(PETSC PETSc>=3.17 IMPORTED_TARGET REQUIRED)


set(HAVE_PETSC ${PETSC_FOUND})

if(${PETSC_FOUND})
  # By default, the scope of imported targets is the current directory. Below we make it known globally.
  set_target_properties(PkgConfig::PETSC PROPERTIES IMPORTED_GLOBAL TRUE)
  dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_PETSC=1"
                              INCLUDE_DIRS "${PETSC_INCLUDE_DIRS}"
                              LIBRARIES PkgConfig::PETSC
                              )
endif()
