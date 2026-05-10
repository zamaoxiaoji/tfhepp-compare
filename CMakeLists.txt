cmake_minimum_required(VERSION 3.16)
project(tfhe++ CXX C)
set(CMAKE_CXX_STANDARD 20)
option(DEBUG "Debug mode" OFF)
if(DEBUG)
  set(CMAKE_CXX_FLAGS
      "-march=native -g -Wall -Wextra -pedantic -Wno-sign-compare")
else()
  set(CMAKE_CXX_FLAGS
      "-march=native -O3 -g -funroll-loops -Wall -Wextra -pedantic -Wno-sign-compare"
  )
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

option(USE_80BIT_SECURITY "Use 80bit security parameter(faster)" OFF)
option(USE_CGGI19 "Use the parameter set proposed in CGGI19" OFF)
option(USE_CONCRETE "Use the parameter set proposed in CONCRETE" OFF)
option(USE_TFHE_RS "Use the parameter set proposed in TFHE-rs" OFF)
option(USE_TERNARY "Use ternary secret keys" OFF)
option(USE_COMPRESS "Use compressed ciphertexts" OFF)
option(USE_SUBSET_KEY "Use subset key" OFF)
option(USE_HE3DB_COMPAT "Use HE3DB-compatible parameters (lvl0 uint32_t, n=672)" OFF)

option(USE_TERNARY_CMUX "Use ternary cmux" OFF)
option(USE_KEY_BUNDLE "Use key bundle algorithm" OFF)

option(USE_BLAKE3 "Use blake3 as CSPRNG" ON)
option(USE_RANDEN "Use randen as CSPRNG" OFF)

option(USE_AVX512 "Use AVX512 ver. SPQLIOS" OFF)
option(USE_CONCRETE_FFT "Use concrete-fft" OFF)
option(USE_MKL "Use Intel MKL" OFF)
option(USE_FFTW3 "Use FFTW3" OFF)
option(USE_SPQLIOX_AARCH64 "Use spqliox_aarch64" OFF)
option(USE_HEXL "Use Intel HEXL" OFF)

option(ENABLE_TEST "Build tests" OFF)
option(ENABLE_BENCHMARK "Build benchmarks" OFF)
option(ENABLE_TUTORIAL "Build tutorial" OFF)
option(USE_PERF "Use Google Profiler" OFF)

option(ENABLE_SHARED "Build as shared libraries" OFF)

set(TFHEpp_DEFINITIONS
    ""
    PARENT_SCOPE)

if(USE_BLAKE3)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_BLAKE3"
      PARENT_SCOPE)
  add_compile_definitions(USE_BLAKE3)
  add_subdirectory(thirdparties/BLAKE3/c)
elseif(USE_RANDEN)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_RANDEN"
      PARENT_SCOPE)
  add_compile_definitions(USE_RANDEN)
endif()

if(USE_HE3DB_COMPAT)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_HE3DB_COMPAT"
      PARENT_SCOPE)
  add_compile_definitions(USE_HE3DB_COMPAT)
elseif(USE_80BIT_SECURITY)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_80BIT_SECURITY"
      PARENT_SCOPE)
  add_compile_definitions(USE_80BIT_SECURITY)
elseif(USE_COMPRESS)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_COMPRESS"
      PARENT_SCOPE)
  add_compile_definitions(USE_COMPRESS)
elseif(USE_CGGI19)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_CGGI19"
      PARENT_SCOPE)
  add_compile_definitions(USE_CGGI19)
elseif(USE_CONCRETE)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_CONCRETE"
      PARENT_SCOPE)
  add_compile_definitions(USE_CONCRETE)
elseif(USE_TFHE_RS)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_TFHE_RS"
      PARENT_SCOPE)
  add_compile_definitions(USE_TFHE_RS)
elseif(USE_TERNARY_CMUX)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_TERNARY_CMUX;USE_TERNARY"
      PARENT_SCOPE)
  add_compile_definitions(USE_TERNARY)
  add_compile_definitions(USE_TERNARY_CMUX)
elseif(USE_TERNARY)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_TERNARY"
      PARENT_SCOPE)
  add_compile_definitions(USE_TERNARY)
endif()

if(NOT USE_TERNARY)
  if(USE_KEY_BUNDLE)
    set(TFHEpp_DEFINITIONS
        "${TFHEpp_DEFINITIONS};USE_KEY_BUNDLE"
        PARENT_SCOPE)
    add_compile_definitions(USE_KEY_BUNDLE)
  endif()
endif()

# if(USE_AVX512) string(APPEND CMAKE_CXX_FLAGS " -mprefer-vector-width=512")
# endif()

if(USE_FFTW3)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_FFTW3"
      PARENT_SCOPE)
  add_compile_definitions(USE_FFTW3)
  add_subdirectory(thirdparties/fftw)
elseif(USE_MKL)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_MKL;USE_INTERLEAVED_FORMAT"
      PARENT_SCOPE)
  add_compile_definitions(USE_MKL)
  add_compile_definitions(USE_INTERLEAVED_FORMAT)
  find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
  include_directories(${MKLROOT}/include)
  add_subdirectory(thirdparties/mkl)
elseif(USE_SPQLIOX_AARCH64)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_SPQLIOX_AARCH64"
      PARENT_SCOPE)
  add_compile_definitions(USE_SPQLIOX_AARCH64)
  add_subdirectory(thirdparties/spqliox_aarch64)
  # Check if the platform is macOS and the architecture is ARM64
  if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
    include_directories(/opt/homebrew/include)
  endif()
elseif(USE_CONCRETE_FFT)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_CONCRETE_FFT;USE_INTERLEAVED_FORMAT"
      PARENT_SCOPE)
  add_compile_definitions(USE_CONCRETE_FFT)
  add_compile_definitions(USE_INTERLEAVED_FORMAT)
  add_subdirectory(thirdparties/concrete-fft)
else()
  add_subdirectory(thirdparties/spqlios)
endif()

if(USE_HEXL)
  set(TFHEpp_DEFINITIONS
      "${TFHEpp_DEFINITIONS};USE_HEXL"
      PARENT_SCOPE)
  add_compile_definitions(USE_HEXL)
  add_subdirectory(thirdparties/hexl)
  # set(CMAKE_CXX_FLAGS "-march=native -O3 -g -funroll-loops -Wall -Wextra
  # -pedantic -Wno-sign-compare -mprefer-vector-width=512" )
endif()

# For OpenMP
find_package(OpenMP)
if(OpenMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if(USE_RANDEN)
  add_subdirectory(thirdparties/randen)
endif()
add_subdirectory(src)

function(add_my_ethmsb_target target source)
  add_executable(${target} EXCLUDE_FROM_ALL ${source})
  target_include_directories(${target} PRIVATE ${PROJECT_SOURCE_DIR})
  target_link_libraries(${target} tfhe++)
  target_compile_features(${target} PRIVATE cxx_std_20)
endfunction()

set(HE3DB_SOURCE_DIR "${PROJECT_SOURCE_DIR}/../HE3DB/src")
set(MY_ETHMSB_HE3DB_SOURCES
    "${HE3DB_SOURCE_DIR}/HEDB/comparison/tfhepp_utils.cpp"
    "${HE3DB_SOURCE_DIR}/HEDB/comparison/operators.cpp"
    "${HE3DB_SOURCE_DIR}/HEDB/comparison/extract_msb.cpp")

function(add_my_ethmsb_he3db_target target source)
  add_executable(${target} EXCLUDE_FROM_ALL ${source}
                                             ${MY_ETHMSB_HE3DB_SOURCES})
  target_include_directories(${target} PRIVATE ${PROJECT_SOURCE_DIR}
                                               ${HE3DB_SOURCE_DIR}
                                               ${HE3DB_SOURCE_DIR}/HEDB/comparison)
  target_link_libraries(${target} tfhe++)
  target_compile_features(${target} PRIVATE cxx_std_20)
endfunction()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_tests.cpp")
  add_my_ethmsb_target(my_ethmsb_tests my_ethmsb_tests.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_bench.cpp")
  add_my_ethmsb_target(my_ethmsb_bench my_ethmsb_bench.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_static_tests.cpp")
  add_my_ethmsb_target(my_ethmsb_static_tests my_ethmsb_static_tests.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_trace_one.cpp")
  add_my_ethmsb_target(my_ethmsb_trace_one my_ethmsb_trace_one.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_pbs_matrix_tests.cpp")
  add_my_ethmsb_target(my_ethmsb_pbs_matrix_tests
                       my_ethmsb_pbs_matrix_tests.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_switchband.cpp")
  add_my_ethmsb_target(my_ethmsb_switchband my_ethmsb_switchband.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_margin_cert.cpp")
  add_my_ethmsb_target(my_ethmsb_margin_cert my_ethmsb_margin_cert.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_random_resumable.cpp")
  add_my_ethmsb_target(my_ethmsb_random_resumable
                       my_ethmsb_random_resumable.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_iks_precision_tests.cpp")
  add_my_ethmsb_target(my_ethmsb_iks_precision_tests
                       my_ethmsb_iks_precision_tests.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_h3compat_targeted_tests.cpp")
  add_my_ethmsb_target(my_ethmsb_h3compat_targeted_tests
                       my_ethmsb_h3compat_targeted_tests.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_he3db_fairness_audit.cpp")
  add_my_ethmsb_target(my_ethmsb_he3db_fairness_audit
                       my_ethmsb_he3db_fairness_audit.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_offset_semantics_tests.cpp")
  add_my_ethmsb_target(my_ethmsb_offset_semantics_tests
                       my_ethmsb_offset_semantics_tests.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_delta_prime_semantics_tests.cpp")
  add_my_ethmsb_target(my_ethmsb_delta_prime_semantics_tests
                       my_ethmsb_delta_prime_semantics_tests.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_fair_he3db_comparison_test.cpp")
  add_my_ethmsb_target(my_ethmsb_fair_he3db_comparison_test
                       my_ethmsb_fair_he3db_comparison_test.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_he3db_style_tests.cpp")
  add_my_ethmsb_target(my_ethmsb_he3db_style_tests
                       my_ethmsb_he3db_style_tests.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_pruned_gapmsb_tests.cpp")
  add_my_ethmsb_target(my_ethmsb_pruned_gapmsb_tests
                       my_ethmsb_pruned_gapmsb_tests.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_pruned_cmux_audit.cpp")
  add_my_ethmsb_target(my_ethmsb_pruned_cmux_audit
                       my_ethmsb_pruned_cmux_audit.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_pruned_margin_cert.cpp")
  add_my_ethmsb_target(my_ethmsb_pruned_margin_cert
                       my_ethmsb_pruned_margin_cert.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_pruned_failure_decomp.cpp")
  add_my_ethmsb_target(my_ethmsb_pruned_failure_decomp
                       my_ethmsb_pruned_failure_decomp.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_pruned_failure_sweep.cpp")
  add_my_ethmsb_target(my_ethmsb_pruned_failure_sweep
                       my_ethmsb_pruned_failure_sweep.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_pruned_gapmsb_bench_smoke.cpp")
  add_my_ethmsb_target(my_ethmsb_pruned_gapmsb_bench_smoke
                       my_ethmsb_pruned_gapmsb_bench_smoke.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_pruned_tex_claim_audit.cpp")
  add_my_ethmsb_target(my_ethmsb_pruned_tex_claim_audit
                       my_ethmsb_pruned_tex_claim_audit.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_direct_mask_periodicity_audit.cpp")
  add_my_ethmsb_target(my_ethmsb_direct_mask_periodicity_audit
                       my_ethmsb_direct_mask_periodicity_audit.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_fused_lut_truth_table_audit.cpp")
  add_my_ethmsb_target(my_ethmsb_fused_lut_truth_table_audit
                       my_ethmsb_fused_lut_truth_table_audit.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_fused_final_tests.cpp")
  add_my_ethmsb_target(my_ethmsb_fused_final_tests
                       my_ethmsb_fused_final_tests.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_fused_vs_conversion_ab.cpp")
  add_my_ethmsb_target(my_ethmsb_fused_vs_conversion_ab
                       my_ethmsb_fused_vs_conversion_ab.cpp)
endif()

if(EXISTS
   "${PROJECT_SOURCE_DIR}/my_ethmsb_pruned_fused_final_bench_smoke.cpp")
  add_my_ethmsb_target(my_ethmsb_pruned_fused_final_bench_smoke
                       my_ethmsb_pruned_fused_final_bench_smoke.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_compare_he3db.cpp"
   AND EXISTS "${HE3DB_SOURCE_DIR}/HEDB/comparison/HomCompare.h")
  add_my_ethmsb_he3db_target(my_ethmsb_compare_he3db
                             my_ethmsb_compare_he3db.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_bench_he3db_style.cpp"
   AND EXISTS "${HE3DB_SOURCE_DIR}/HEDB/comparison/HomCompare.h")
  add_my_ethmsb_he3db_target(my_ethmsb_bench_he3db_style
                             my_ethmsb_bench_he3db_style.cpp)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_vs_he3db_test.cpp")
  add_executable(my_ethmsb_vs_he3db_test EXCLUDE_FROM_ALL
      my_ethmsb_vs_he3db_test.cpp
      ethmsb.cpp)
  target_include_directories(my_ethmsb_vs_he3db_test PRIVATE
      ${PROJECT_SOURCE_DIR})
  target_link_libraries(my_ethmsb_vs_he3db_test tfhe++)
  target_compile_features(my_ethmsb_vs_he3db_test PRIVATE cxx_std_20)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/my_ethmsb_pruned_fast_test.cpp")
  add_executable(my_ethmsb_pruned_fast_test EXCLUDE_FROM_ALL
      my_ethmsb_pruned_fast_test.cpp
      ethmsb.cpp)
  target_include_directories(my_ethmsb_pruned_fast_test PRIVATE
      ${PROJECT_SOURCE_DIR})
  target_link_libraries(my_ethmsb_pruned_fast_test tfhe++)
  target_compile_features(my_ethmsb_pruned_fast_test PRIVATE cxx_std_20)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/experimental/weighted_bitapprox_probe.cpp")
  add_executable(my_weighted_bitapprox_probe EXCLUDE_FROM_ALL
      experimental/weighted_bitapprox_probe.cpp)
  target_include_directories(my_weighted_bitapprox_probe PRIVATE
      ${PROJECT_SOURCE_DIR}
      ${PROJECT_SOURCE_DIR}/include)
  target_link_libraries(my_weighted_bitapprox_probe tfhe++)
  target_compile_features(my_weighted_bitapprox_probe PRIVATE cxx_std_20)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/experimental/weighted_bitapprox_probe.cpp")
  add_executable(my_pruned_weightedapprox_full_pipeline_probe EXCLUDE_FROM_ALL
      experimental/weighted_bitapprox_probe.cpp)
  target_include_directories(my_pruned_weightedapprox_full_pipeline_probe PRIVATE
      ${PROJECT_SOURCE_DIR}
      ${PROJECT_SOURCE_DIR}/include)
  target_link_libraries(my_pruned_weightedapprox_full_pipeline_probe tfhe++)
  target_compile_features(my_pruned_weightedapprox_full_pipeline_probe PRIVATE
      cxx_std_20)
endif()

if(ENABLE_TEST)
  add_subdirectory(test)
  add_subdirectory(thirdparties/AES)
endif()
if(ENABLE_BENCHMARK)
  add_subdirectory(benchmark)
endif()
if(ENABLE_TUTORIAL)
  add_subdirectory(tutorial)
endif()

install(TARGETS tfhe++ LIBRARY DESTINATION lib)
