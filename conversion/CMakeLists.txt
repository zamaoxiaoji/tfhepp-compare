add_library(tfhepp_ckks_repack STATIC
    ckks_repack.cpp
    ckks_relational.cpp
    ckks_rns_masks.cpp)

target_include_directories(tfhepp_ckks_repack PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(tfhepp_ckks_repack PUBLIC tfhe++ seal)
target_compile_features(tfhepp_ckks_repack PUBLIC cxx_std_20)
