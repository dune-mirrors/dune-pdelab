#=========================
# Linear Matrix-Free Tests
#=========================

dune_add_test(NAME matrix_free_linear_partial_matrix_free_fastdg
              SOURCES matrix_free_linear.cc
              COMPILE_DEFINITIONS PARTIAL_MATRIX_FREE FASTDG
              CMAKE_GUARD HAVE_EIGEN)
add_dune_eigen_flags(matrix_free_linear_partial_matrix_free_fastdg)

dune_add_test(NAME matrix_free_linear_partial_matrix_free_nonfast
              SOURCES matrix_free_linear.cc
              COMPILE_DEFINITIONS PARTIAL_MATRIX_FREE
              CMAKE_GUARD HAVE_EIGEN)
add_dune_eigen_flags(matrix_free_linear_partial_matrix_free_nonfast)

dune_add_test(NAME matrix_free_linear_fully_matrix_free_fastdg
              SOURCES matrix_free_linear.cc
              COMPILE_DEFINITIONS FULLY_MATRIX_FREE FASTDG
              CMAKE_GUARD HAVE_EIGEN)
add_dune_eigen_flags(matrix_free_linear_fully_matrix_free_fastdg)

dune_add_test(NAME matrix_free_linear_fully_matrix_free_nonfast
              SOURCES matrix_free_linear.cc
              COMPILE_DEFINITIONS FULLY_MATRIX_FREE
              CMAKE_GUARD HAVE_EIGEN)
add_dune_eigen_flags(matrix_free_linear_fully_matrix_free_nonfast)

dune_add_test(NAME matrix_free_linear_matrix_free_sor_fastdg
              SOURCES matrix_free_linear.cc
              COMPILE_DEFINITIONS SOR FASTDG
              CMAKE_GUARD HAVE_EIGEN)
add_dune_eigen_flags(matrix_free_linear_matrix_free_sor_fastdg)

# SOR only works with FASTDG
# dune_add_test(NAME matrix_free_linear_matrix_free_sor_nonfast
#               SOURCES matrix_free_linear.cc
#               COMPILE_DEFINITIONS SOR)

dune_symlink_to_source_files(FILES matrix_free_linear.ini)

#============================
# Nonlinear Matrix-Free Tests
#============================
dune_add_test(NAME matrix_free_nonlinear_partial_matrix_free_fastdg
              SOURCES matrix_free_nonlinear.cc
              COMPILE_DEFINITIONS PARTIAL_MATRIX_FREE FASTDG
              CMAKE_GUARD HAVE_EIGEN)
add_dune_eigen_flags(matrix_free_nonlinear_partial_matrix_free_fastdg)

dune_add_test(NAME matrix_free_nonlinear_partial_matrix_free_nonfast
              SOURCES matrix_free_nonlinear.cc
              COMPILE_DEFINITIONS PARTIAL_MATRIX_FREE
              CMAKE_GUARD HAVE_EIGEN)
add_dune_eigen_flags(matrix_free_nonlinear_partial_matrix_free_nonfast)

dune_add_test(NAME matrix_free_nonlinear_fully_matrix_free_fastdg
              SOURCES matrix_free_nonlinear.cc
              COMPILE_DEFINITIONS FULLY_MATRIX_FREE FASTDG
              CMAKE_GUARD HAVE_EIGEN)
add_dune_eigen_flags(matrix_free_nonlinear_fully_matrix_free_fastdg)

dune_add_test(NAME matrix_free_nonlinear_fully_matrix_free_nonfast
              SOURCES matrix_free_nonlinear.cc
              COMPILE_DEFINITIONS FULLY_MATRIX_FREE
              CMAKE_GUARD HAVE_EIGEN)
add_dune_eigen_flags(matrix_free_nonlinear_fully_matrix_free_nonfast)

dune_add_test(NAME matrix_free_nonlinear_matrix_free_sor_fastdg
              SOURCES matrix_free_nonlinear.cc
              COMPILE_DEFINITIONS SOR FASTDG
              CMAKE_GUARD HAVE_EIGEN)
add_dune_eigen_flags(matrix_free_nonlinear_matrix_free_sor_fastdg)

# SOR only works with FASTDG
# dune_add_test(NAME matrix_free_nonlinear_matrix_free_sor_nonfast
#               SOURCES matrix_free_nonlinear.cc
#               COMPILE_DEFINITIONS SOR)

dune_symlink_to_source_files(FILES matrix_free_nonlinear.ini)
