# Homomorphic comparison library:
#   - Algorithm 1: ETHMSB + gap offset                    (samplepaper.tex)
#   - Algorithm 2: Pruned BitExtract + B2A + ETHMSB       (Chapter 3)
add_library(tfhepp_compare STATIC
    tfhepp_utils.cpp
    ethmsb.cpp
    bitextract_qhalf.cpp
    pruned_three_pbs.cpp
    micro_pbs.cpp)

target_include_directories(tfhepp_compare PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(tfhepp_compare PUBLIC tfhe++)
target_compile_features(tfhepp_compare PUBLIC cxx_std_20)
