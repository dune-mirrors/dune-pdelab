dune_add_test(SOURCES testtopologicalinterleaving.cc)
dune_add_test(SOURCES testlagrangebasis.cc)

# Test for comparing testlagrangebasis output files to reference files
set(TEST_LAGRANGE_INDEX_TREE_OUTPUT_FILES
    blockedByEntity_Q2_index_tree.txt
    composite_power5_flatByEntity_power2_blockedByEntity_flatLexicographic_Q2_index_tree.txt
    composite_power5_flatByEntity_power2_blockedByEntity_Q2_index_tree.txt
    composite_power5_lagrange2_blockedLexicographic_power2_blockedByEntity_blockedLexicographic_Q2_index_tree.txt
    composite_power5_lagrange2_flatLexicographic_power2_flatByEntity_flatLexicographic_Q2_index_tree.txt
    composite_power5_power2_blockedByEntity_flatByEntity_Q2_index_tree.txt
    composite_power5_power2_blockedByEntity_Q1_index_tree.txt
    composite_power5_power2_flatByEntity_Q2_index_tree.txt
    flatByEntity_Q2_index_tree.txt
    power2_blockedByEntity_blockedLexicographic_Q2_index_tree.txt
    power2_blockedByEntity_Q2_index_tree.txt
    power5_blockedByEntity_flatByEntity_Q2_index_tree.txt
    power5_blockedByEntity_flatLexicographic_Q2_index_tree.txt
    power5_flatByEntity_blockedByEntity_Q2_index_tree.txt
    power5_flatByEntity_blockedLexicographic_Q2_index_tree.txt
    power5_flatByEntity_flatLexicographic_Q2_index_tree.txt
    power5_flatByEntity_Q2_index_tree.txt
    TaylorHood_blockedByEntity_blockedByEntity_index_tree.txt
    TaylorHood_blockedByEntity_flatByEntity_index_tree.txt
    TaylorHood_blockedLexicographic_flatByEntity_index_tree.txt
    TaylorHood_flatByEntity_index_tree.txt
)

foreach(output_file IN LISTS TEST_LAGRANGE_INDEX_TREE_OUTPUT_FILES)
    add_test(
        NAME compare_${output_file}
        COMMAND ${CMAKE_COMMAND} -E compare_files
            ${CMAKE_CURRENT_BINARY_DIR}/${output_file}
            ${CMAKE_CURRENT_SOURCE_DIR}/reference/${output_file}
    )
    set_tests_properties(compare_${output_file} PROPERTIES DEPENDS testlagrangebasis)
endforeach()
