#include "algorithms.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

std::vector<double> DecryptSlots(
    const PaperReview::CkksRuntime& runtime,
    const PaperReview::CkksCiphertext& ct,
    std::size_t slots) {
    lbcrypto::Plaintext plain;
    runtime.cc->Decrypt(runtime.keys.secretKey, ct, &plain);
    plain->SetLength(slots);
    return plain->GetRealPackedValue();
}

void RequireClose(
    const std::vector<double>& actual,
    const std::vector<double>& expected,
    double tolerance,
    const char* label) {
    if (actual.size() < expected.size())
        throw std::runtime_error(std::string(label) + ": not enough slots");
    for (std::size_t i = 0; i < expected.size(); i++) {
        if (std::abs(actual[i] - expected[i]) > tolerance) {
            throw std::runtime_error(std::string(label) + ": slot mismatch");
        }
    }
}

}  // namespace

int main() {
    constexpr std::size_t slots = 4;
    auto runtime = PaperReview::MakeCkksRuntime(slots, 6);

    const std::vector<PaperReview::AttributeSpec> attrs{{2}};
    const auto domain = PaperReview::EnumerateDomainTuples(attrs);
    const auto basis = PaperReview::EnumerateMonomialBasis(attrs);
    const auto matrix = PaperReview::BuildBasisMatrix(domain, basis);
    const auto alpha = PaperReview::SolveCoeffTableForAllTargets(matrix);
    if (!PaperReview::VerifyAlphaTable(matrix, alpha))
        throw std::runtime_error("AlphaTable verification failed");

    const auto group_key = PaperReview::EncryptVector(
        runtime.cc, runtime.keys.publicKey, {0.0, 1.0, 0.0, 1.0}, slots);
    const auto value_col = PaperReview::EncryptVector(
        runtime.cc, runtime.keys.publicKey, {10.0, 20.0, 30.0, 40.0}, slots);
    const auto grouped = PaperReview::MatrixGroupBySum(
        runtime.cc, runtime.keys.publicKey, attrs, {group_key},
        value_col, basis, alpha, slots);
    RequireClose(DecryptSlots(runtime, grouped[0], slots), {40.0, 40.0, 40.0, 40.0}, 0.25, "GROUP BY key 0");
    RequireClose(DecryptSlots(runtime, grouped[1], slots), {60.0, 60.0, 60.0, 60.0}, 0.25, "GROUP BY key 1");

    const auto left_key = PaperReview::EncryptVector(
        runtime.cc, runtime.keys.publicKey, {0.0, 1.0, 0.0, 1.0}, slots);
    const auto right_key = PaperReview::EncryptVector(
        runtime.cc, runtime.keys.publicKey, {0.0, 1.0, 0.0, 0.0}, slots);
    const auto right_payload = PaperReview::EncryptVector(
        runtime.cc, runtime.keys.publicKey, {7.0, 9.0, 0.0, 0.0}, slots);
    const auto joined = PaperReview::LookupJoin(
        runtime.cc, runtime.keys.publicKey, left_key, {},
        right_key, {right_payload}, basis, alpha, slots);
    RequireClose(DecryptSlots(runtime, joined[0], slots), {7.0, 9.0, 7.0, 9.0}, 0.5, "LOOKUP JOIN");

    std::cout << "CKKS GROUP BY/JOIN smoke: ok" << std::endl;
    return 0;
}
