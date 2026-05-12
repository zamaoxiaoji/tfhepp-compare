#include <iostream>
#include <vector>

#include "ckks_repack.h"
#include "params.hpp"

int main()
{
    using P = TFHEpp::lvl1param;

    auto parms = tfhepp_ckks::MakeDefaultCKKSParameters();
    auto config = tfhepp_ckks::DefaultRepackConfig<P>();

    using GenFn = void (*)(tfhepp_ckks::RepackEvaluationKey &,
                           const TFHEpp::SecretKey &, double,
                           seal::CKKSEncoder &, const seal::Encryptor &,
                           const seal::SEALContext &);
    using PackFn = void (*)(seal::Ciphertext &,
                            const std::vector<TFHEpp::TLWE<P>> &,
                            const tfhepp_ckks::RepackEvaluationKey &,
                            const tfhepp_ckks::RepackConfig &,
                            seal::CKKSEncoder &, const seal::GaloisKeys &,
                            seal::RelinKeys &, seal::Evaluator &,
                            seal::SEALContext &);

    GenFn gen = &tfhepp_ckks::GenerateRepackKey<P>;
    PackFn pack = &tfhepp_ckks::PackLWEsToCKKS<P>;

    if (parms.poly_modulus_degree() != 65536 || config.message_scale <= 0.0 ||
        config.lwe_modulus <= 0.0 || config.key_scale <= 0.0 || gen == nullptr ||
        pack == nullptr) {
        return 1;
    }

    std::cout << "ckks_repack_interface_smoke ok" << std::endl;
    return 0;
}
