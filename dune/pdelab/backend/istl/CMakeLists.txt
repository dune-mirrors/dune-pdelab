add_subdirectory(geneo)
add_subdirectory(matrixfree)

install(FILES
  bcrsmatrix.hh
  bcrsmatrixbackend.hh
  bcrspattern.hh
  blockmatrixdiagonal.hh
  cg_to_dg_prolongation.hh
  descriptors.hh
  dunefunctions.hh
  forwarddeclarations.hh
  istlsolverbackend.hh
  matrixhelpers.hh
  novlpistlsolverbackend.hh
  ovlp_amg_dg_backend.hh
  ovlpistlsolverbackend.hh
  parallelhelper.hh
  patternstatistics.hh
  seq_amg_dg_backend.hh
  seqistlsolverbackend.hh
  tags.hh
  utility.hh
  vector.hh
  vectorhelpers.hh
  vectoriterator.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/backend/istl)
