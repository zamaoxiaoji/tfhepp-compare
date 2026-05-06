find_package(OpenFHE CONFIG REQUIRED)
find_package(OpenMP)

add_library(paper_native_baselines native_baseline_experiments.cpp)
target_include_directories(paper_native_baselines PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(paper_native_baselines PUBLIC tfhe++)
target_compile_features(paper_native_baselines PUBLIC cxx_std_20)
if(OpenMP_CXX_FOUND)
    target_link_libraries(paper_native_baselines PUBLIC OpenMP::OpenMP_CXX)
endif()

add_library(paper_algorithms algorithms.cpp tpch_experiments.cpp chap3_compare_experiment.cpp)
target_include_directories(paper_algorithms PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/include
)
target_include_directories(paper_algorithms SYSTEM PUBLIC
    ${OpenFHE_INCLUDE}
    ${OpenFHE_INCLUDE}/third-party/include
    ${OpenFHE_INCLUDE}/core
    ${OpenFHE_INCLUDE}/pke
    ${OpenFHE_INCLUDE}/binfhe
)
target_link_directories(paper_algorithms PUBLIC ${OpenFHE_LIBDIR})
target_link_libraries(paper_algorithms PUBLIC tfhe++ ${OpenFHE_SHARED_LIBRARIES})
target_compile_features(paper_algorithms PUBLIC cxx_std_20)
if(OpenMP_CXX_FOUND)
    target_link_libraries(paper_algorithms PUBLIC OpenMP::OpenMP_CXX)
endif()

if(UNIX AND NOT APPLE)
    target_link_options(paper_algorithms PUBLIC -Wl,--no-as-needed)
endif()

macro(add_paper_experiment name source)
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE paper_algorithms paper_native_baselines)
endmacro()

macro(add_chap5_reference_experiment name source)
    add_executable(${name} ${source})
    target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_compile_features(${name} PRIVATE cxx_std_17)
endmacro()

add_paper_experiment(exp_chap3_rel_comp exp_chap3_rel_comp.cpp)
add_paper_experiment(exp_chap3_accuracy exp_chap3_accuracy.cpp)
add_paper_experiment(exp_chap3_boundary exp_chap3_boundary.cpp)
add_paper_experiment(exp_chap3_k_impact exp_chap3_k_impact.cpp)
add_paper_experiment(exp_chap3_compare_test exp_chap3_compare_test.cpp)

add_paper_experiment(exp_chap4_group_by exp_chap4_group_by.cpp)
add_paper_experiment(exp_chap4_join exp_chap4_join.cpp)
add_paper_experiment(exp_chap4_rns exp_chap4_rns.cpp)

add_paper_experiment(exp_chap5_query6 exp_chap5_query6.cpp)
add_paper_experiment(exp_chap5_query14 exp_chap5_query14.cpp)
add_paper_experiment(exp_chap5_query3 exp_chap5_query3.cpp)
add_paper_experiment(exp_chap5_query5 exp_chap5_query5.cpp)

add_paper_experiment(review_repack_compile_check review_repack_compile_check.cpp)
add_paper_experiment(review_ckks_operator_smoke review_ckks_operator_smoke.cpp)
