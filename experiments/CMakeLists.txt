macro(add_ours_experiment name source)
    add_executable(${name} ${source})
    target_include_directories(${name} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${PROJECT_SOURCE_DIR}/paper_experiments
        ${PROJECT_SOURCE_DIR}/include
    )
    target_link_libraries(${name} PRIVATE paper_algorithms)
    target_compile_features(${name} PRIVATE cxx_std_20)
endmacro()

set(HE3DB_ROOT "/home/wrn/HE3DB" CACHE PATH "Path to a local HE3DB checkout")

macro(add_he3db_latency_experiment name source)
    if(EXISTS "${HE3DB_ROOT}/src/HEDB/comparison/comparison.h"
       AND EXISTS "${HE3DB_ROOT}/build/lib/libHEDB.a")
        add_executable(${name} ${source})
        target_include_directories(${name} BEFORE PRIVATE
            ${HE3DB_ROOT}/src
            ${HE3DB_ROOT}/thirdparty/TFHEpp/include
            ${HE3DB_ROOT}/thirdparty/TFHEpp/thirdparties/cereal/include
            ${HE3DB_ROOT}/thirdparty/TFHEpp/thirdparties/randen
            ${HE3DB_ROOT}/thirdparty/TFHEpp/thirdparties/spqlios
        )
        target_compile_definitions(${name} PRIVATE USE_RANDEN)
        target_compile_features(${name} PRIVATE cxx_std_20)
        target_link_libraries(${name} PRIVATE
            ${HE3DB_ROOT}/build/lib/libHEDB.a
            ${HE3DB_ROOT}/build/lib/libtfhe++.a
            ${HE3DB_ROOT}/build/lib/libspqlios.a
            ${HE3DB_ROOT}/build/lib/libranden.a
        )
        if(OpenMP_CXX_FOUND)
            target_link_libraries(${name} PRIVATE OpenMP::OpenMP_CXX)
        endif()
    else()
        message(STATUS "Skipping ${name}: HE3DB_ROOT=${HE3DB_ROOT} is missing headers or build/lib/libHEDB.a")
    endif()
endmacro()

add_ours_experiment(exp_correctness exp_correctness.cpp)
add_ours_experiment(exp_boundary exp_boundary.cpp)
add_ours_experiment(exp_k_choice exp_k_choice.cpp)
add_ours_experiment(exp_latency_relation exp_latency_relation.cpp)
add_ours_experiment(exp_latency_ops exp_latency_ops.cpp)
add_ours_experiment(exp_pruning_counter exp_pruning_counter.cpp)
add_ours_experiment(exp_recursive_delta_check exp_recursive_delta_check.cpp)
add_ours_experiment(exp_recursive_schedule_check exp_recursive_schedule_check.cpp)
add_ours_experiment(exp_recursive_gap32_check exp_recursive_gap32_check.cpp)
add_he3db_latency_experiment(exp_he3db_compare_latency exp_he3db_compare_latency.cpp)
