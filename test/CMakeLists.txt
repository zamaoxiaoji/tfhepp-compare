file(GLOB test_sources RELATIVE "${CMAKE_CURRENT_LIST_DIR}" "*.cpp")
# comparison_test needs to link tfhepp_compare in addition to tfhe++/AES,
# so it is built explicitly below rather than through the auto-glob loop.
list(REMOVE_ITEM test_sources "comparison_test.cpp")

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

# comparison/ benchmark — exercises both ETHMSB and Pruned-3-PBS.
if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/comparison_test.cpp" AND TARGET tfhepp_compare)
    add_executable(comparison_test comparison_test.cpp)
    target_link_libraries(comparison_test tfhepp_compare tfhe++ AES)
    target_compile_features(comparison_test PRIVATE cxx_std_20)
endif()

# For GMP used in INTorus test
find_library(gmp NAMES gmp gmpxx)
if(gmp)
    add_compile_definitions(ENABLE_LSHTEST)
    target_link_libraries(INTorus gmp)
endif()
