add_subdirectory(functions)

install(FILES constraints.hh
              constraintsparameters.hh
              constraintstransformation.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/constraints/common)
