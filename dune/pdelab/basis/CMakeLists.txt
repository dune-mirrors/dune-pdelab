add_subdirectory(protobasis)
add_subdirectory(topologicassociativityforest)

install(FILES merging_strategy.hh
              protobasis.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/basis)
