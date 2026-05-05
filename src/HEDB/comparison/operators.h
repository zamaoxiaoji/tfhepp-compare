#pragma once
#include "cloudkey.hpp"
#include "HEDB/utils/types.h"

namespace HEDB
{
    // HomAND: c0 + c1 - Q/8, then bootstrap
    void HomAND(TLWELvl1 &res, const TLWELvl1 &ca, const TLWELvl1 &cb,
                const TFHEEvalKey &ek, bool result_type);

    // HomOR: c0 + c1 + Q/8, then bootstrap
    void HomOR(TLWELvl1 &res, const TLWELvl1 &ca, const TLWELvl1 &cb,
               const TFHEEvalKey &ek, bool result_type);

    // HomNOT: negate all elements
    template <typename P>
    inline void HomNOT(TFHEpp::TLWE<P> &res, const TFHEpp::TLWE<P> &tlwe)
    {
        for (int i = 0; i <= P::k * P::n; i++) res[i] = -tlwe[i];
    }
} // namespace HEDB
