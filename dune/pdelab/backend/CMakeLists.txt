add_subdirectory(common)
add_subdirectory(eigen)
add_subdirectory(istl)
add_subdirectory(simple)
add_subdirectory(petsc)

install(FILES petsc.hh
              eigen.hh
              interface.hh
              istl.hh
              simple.hh
              solver.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/backend)
