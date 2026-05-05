#include "algorithms.hpp"

#include <cstdint>
#include <vector>

int main() {
    static_assert(PaperReview::kRequiredSecurityBits == 128);

    using LweP = TFHEpp::lvl1param;

    auto key_converter = &PaperReview::TFHEppKeyToOpenFHE<LweP>;
    auto ciphertext_converter = &PaperReview::TFHEppLWEsToOpenFHE<LweP>;
    auto repack_setup = &PaperReview::RepackSetup<LweP>;
    auto repack_execute = &PaperReview::RepackExecute<LweP>;

    (void)key_converter;
    (void)ciphertext_converter;
    (void)repack_setup;
    (void)repack_execute;

    return 0;
}
