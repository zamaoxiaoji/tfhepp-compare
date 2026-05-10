#pragma once
/**
 * @file comparison.h
 * @brief Umbrella header for the TFHEpp comparison library.
 *
 * Bring in:
 *   - tfhepp_compare::ethmsb     — Algorithm 1: ETHMSB + gap offset
 *                                   (samplepaper.tex)
 *   - tfhepp_compare::three_pbs  — Algorithm 2: periodic-pruned BitExtract +
 *                                   B2A + ETHMSB (Chapter 3 of the thesis)
 *
 * Each namespace exposes:
 *   - HomMSB(res, ct, plain_bits, ek, result_type)
 *   - greater_than / greater_than_equal / less_than / less_than_equal / equal
 *
 * Plus the shared utilities in tfhepp_compare:: (encrypt/decrypt helpers,
 * LOGIC/ARITHMETIC selectors, ARI↔LOG conversion, gate operators).
 */
#include "ethmsb.h"
#include "pruned_three_pbs.h"
#include "HomCompare.h"
#include "tfhepp_utils.h"
