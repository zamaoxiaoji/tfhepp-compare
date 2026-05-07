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

add_ours_experiment(exp_correctness exp_correctness.cpp)
add_ours_experiment(exp_boundary exp_boundary.cpp)
add_ours_experiment(exp_k_choice exp_k_choice.cpp)
add_ours_experiment(exp_latency_relation exp_latency_relation.cpp)
add_ours_experiment(exp_latency_ops exp_latency_ops.cpp)
add_ours_experiment(exp_pruning_counter exp_pruning_counter.cpp)
add_ours_experiment(exp_recursive_delta_check exp_recursive_delta_check.cpp)
add_ours_experiment(exp_recursive_schedule_check exp_recursive_schedule_check.cpp)
