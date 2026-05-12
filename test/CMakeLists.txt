file(GLOB test_sources RELATIVE "${CMAKE_CURRENT_LIST_DIR}" "*.cpp")
# comparison_test and micro_pbs_* tests need to link tfhepp_compare in
# addition to tfhe++/AES, so they are built explicitly below rather than
# through the auto-glob loop.
list(REMOVE_ITEM test_sources "comparison_test.cpp")
list(REMOVE_ITEM test_sources "bitextract_qhalf_boundary_probe.cpp")
list(REMOVE_ITEM test_sources "micro_pbs_isolated_test.cpp")
list(REMOVE_ITEM test_sources "micro_pbs_sweep_test.cpp")
list(REMOVE_ITEM test_sources "micro_pbs_full_pipeline_test.cpp")
list(REMOVE_ITEM test_sources "pruned_three_pbs_targeted_probe.cpp")
list(REMOVE_ITEM test_sources "pruned_three_pbs_kappa_sweep.cpp")
list(REMOVE_ITEM test_sources "pruned_bitclear_lsb_sweep.cpp")
list(REMOVE_ITEM test_sources "ckks_repack_interface_smoke.cpp")

file(COPY test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS OWNER_READ OWNER_EXECUTE)
file(COPY show_info.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS OWNER_READ OWNER_EXECUTE)
file(COPY benchmark_docker.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS OWNER_READ OWNER_EXECUTE)

foreach(test_source ${test_sources})
    string( REPLACE ".cpp" "" test_name ${test_source} )
    add_executable(${test_name} ${test_source})
    if(USE_PERF)
        add_compile_definitions(USE_PERF)
        target_link_libraries(${test_name} tfhe++ AES profiler)
    else()
        target_link_libraries(${test_name} tfhe++ AES)
    endif()
endforeach(test_source ${test_sources})

# comparison/ benchmark — exercises ETHMSB, Pruned-3-PBS, and (optionally)
# the original HE3DB HomMSB pipeline if HE3DB sources are found alongside
# this repo (../HE3DB/src/HEDB/comparison/...).
if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/comparison_test.cpp" AND TARGET tfhepp_compare)
    set(_he3db_dir "${PROJECT_SOURCE_DIR}/../HE3DB/src")
    set(_he3db_extra_sources "")
    set(_he3db_extra_includes "")
    if(EXISTS "${_he3db_dir}/HEDB/comparison/extract_msb.cpp")
        list(APPEND _he3db_extra_sources
            "${_he3db_dir}/HEDB/comparison/tfhepp_utils.cpp"
            "${_he3db_dir}/HEDB/comparison/operators.cpp"
            "${_he3db_dir}/HEDB/comparison/extract_msb.cpp"
            "${_he3db_dir}/HEDB/comparison/HomCompare.cpp")
        list(APPEND _he3db_extra_includes
            "${_he3db_dir}"
            "${_he3db_dir}/HEDB/comparison")
    endif()

    add_executable(comparison_test comparison_test.cpp ${_he3db_extra_sources})
    target_link_libraries(comparison_test tfhepp_compare tfhe++ AES)
    target_compile_features(comparison_test PRIVATE cxx_std_20)
    if(_he3db_extra_sources)
        target_include_directories(comparison_test PRIVATE ${_he3db_extra_includes})
        target_compile_definitions(comparison_test PRIVATE HAVE_HE3DB_ORIGINAL=1)
        message(STATUS "comparison_test: HE3DB original sources found at ${_he3db_dir}")
    else()
        message(STATUS "comparison_test: HE3DB original sources NOT found (skipping HE3DB column)")
    endif()
endif()

if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/bitextract_qhalf_boundary_probe.cpp" AND TARGET tfhepp_compare)
    add_executable(bitextract_qhalf_boundary_probe bitextract_qhalf_boundary_probe.cpp)
    target_link_libraries(bitextract_qhalf_boundary_probe tfhepp_compare tfhe++ AES)
    target_compile_features(bitextract_qhalf_boundary_probe PRIVATE cxx_std_20)
endif()

# Isolated Micro-PBS conversion test (UNSAFE_PERFORMANCE_ONLY).
if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/micro_pbs_isolated_test.cpp" AND TARGET tfhepp_compare)
    add_executable(micro_pbs_isolated_test micro_pbs_isolated_test.cpp)
    target_link_libraries(micro_pbs_isolated_test tfhepp_compare tfhe++ AES)
    target_compile_features(micro_pbs_isolated_test PRIVATE cxx_std_20)
endif()

if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/micro_pbs_sweep_test.cpp" AND TARGET tfhepp_compare)
    add_executable(micro_pbs_sweep_test micro_pbs_sweep_test.cpp)
    target_link_libraries(micro_pbs_sweep_test tfhepp_compare tfhe++ AES)
    target_compile_features(micro_pbs_sweep_test PRIVATE cxx_std_20)
endif()

if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/micro_pbs_full_pipeline_test.cpp" AND TARGET tfhepp_compare)
    add_executable(micro_pbs_full_pipeline_test micro_pbs_full_pipeline_test.cpp)
    target_link_libraries(micro_pbs_full_pipeline_test tfhepp_compare tfhe++ AES)
    target_compile_features(micro_pbs_full_pipeline_test PRIVATE cxx_std_20)
endif()

if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/pruned_three_pbs_targeted_probe.cpp" AND TARGET tfhepp_compare)
    add_executable(pruned_three_pbs_targeted_probe pruned_three_pbs_targeted_probe.cpp)
    target_link_libraries(pruned_three_pbs_targeted_probe tfhepp_compare tfhe++ AES)
    target_compile_features(pruned_three_pbs_targeted_probe PRIVATE cxx_std_20)
endif()

if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/pruned_three_pbs_kappa_sweep.cpp" AND TARGET tfhepp_compare)
    add_executable(pruned_three_pbs_kappa_sweep pruned_three_pbs_kappa_sweep.cpp)
    target_link_libraries(pruned_three_pbs_kappa_sweep tfhepp_compare tfhe++ AES)
    target_compile_features(pruned_three_pbs_kappa_sweep PRIVATE cxx_std_20)
endif()

if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/pruned_bitclear_lsb_sweep.cpp" AND TARGET tfhepp_compare)
    add_executable(pruned_bitclear_lsb_sweep pruned_bitclear_lsb_sweep.cpp)
    target_link_libraries(pruned_bitclear_lsb_sweep tfhepp_compare tfhe++ AES)
    target_compile_features(pruned_bitclear_lsb_sweep PRIVATE cxx_std_20)
endif()

if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/ckks_repack_interface_smoke.cpp" AND TARGET tfhepp_ckks_repack)
    add_executable(ckks_repack_interface_smoke ckks_repack_interface_smoke.cpp)
    target_link_libraries(ckks_repack_interface_smoke tfhepp_ckks_repack tfhe++)
    target_compile_features(ckks_repack_interface_smoke PRIVATE cxx_std_20)
endif()


# For GMP used in INTorus test
find_library(gmp NAMES gmp gmpxx)
if(gmp)
    add_compile_definitions(ENABLE_LSHTEST)
    target_link_libraries(INTorus gmp)
endif()
