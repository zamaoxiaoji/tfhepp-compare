find_package(Eigen3 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
    include(FetchContent)
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
        GIT_SHALLOW TRUE)
    FetchContent_GetProperties(eigen)
    if(NOT eigen_POPULATED)
        FetchContent_Populate(eigen)
    endif()
    add_library(tfhepp_eigen INTERFACE)
    target_include_directories(tfhepp_eigen INTERFACE ${eigen_SOURCE_DIR})
    add_library(Eigen3::Eigen ALIAS tfhepp_eigen)
endif()

add_library(tfhepp_ckks_repack STATIC
    ckks_repack.cpp
    ckks_relational.cpp
    ckks_rns_masks.cpp)

target_include_directories(tfhepp_ckks_repack PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(tfhepp_ckks_repack PUBLIC tfhe++ seal Eigen3::Eigen)
target_compile_features(tfhepp_ckks_repack PUBLIC cxx_std_20)
