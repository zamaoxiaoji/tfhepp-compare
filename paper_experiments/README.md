# Paper Experiment Runbook

This directory is for reviewable paper-facing experiments. The executables do
not print fixed manuscript tables. Each target constructs input data, encrypts
the columns or masks, runs the homomorphic operator under test, decrypts only
the final result, and checks it against a plaintext baseline computed in the
same executable.

## Build

```bash
cmake -S . -B build -DENABLE_TEST=ON
cmake --build build --target paper_algorithms \
  exp_chap3_accuracy exp_chap3_compare_test exp_chap3_rel_comp \
  exp_chap3_boundary exp_chap3_k_impact \
  exp_chap4_group_by exp_chap4_join exp_chap4_rns \
  exp_chap5_query6 exp_chap5_query14 exp_chap5_query3 exp_chap5_query5 \
  review_repack_compile_check review_ckks_operator_smoke
```

## Chapter 3 Comparison Experiments

`exp_chap3_accuracy`, `exp_chap3_compare_test`, `exp_chap3_rel_comp`, and
`exp_chap3_boundary` call `PaperReview::RunChapter3CompareExperiment`, which
performs:

1. TFHEpp secret-key, bootstrapping-key, key-switching-key, and Meta-PBS
   TruncRepeat-key generation.
2. Plaintext case generation for signed-safe unsigned operands. A `bits`
   value means the input operands are in `[0, 2^(bits-1)-1]`; the encrypted
   signed difference uses one extra sign bit.
3. TFHE lvl2 TLWE encryption of `lhs` and `rhs`.
4. `PaperReview::Chapter3HomCompare`, which dispatches to the `metapbs2`
   recursive HomMSB/GapMSB comparison path.
5. lvl1 sign decryption and row-by-row comparison with the plaintext predicate
   `lhs > rhs`.

`exp_chap3_k_impact` follows the Chapter 3 zeroed-bit experiment directly: it
encrypts boundary-neighborhood messages, runs `GapMSB_k` for each configured
zeroed bit `k`, and reports `w_k`, the first-round pruning period, pruning
statistics, latency, and the decrypted MSB check.

Examples:

```bash
./build/paper_experiments/exp_chap3_accuracy
./build/paper_experiments/exp_chap3_compare_test --full
./build/paper_experiments/exp_chap3_boundary --bits 16 --radius 4
./build/paper_experiments/exp_chap3_k_impact --sweep --radius 4
```

## Chapter 4 CKKS Operators

`exp_chap4_group_by`, `exp_chap4_join`, and `exp_chap4_rns` are operator
experiments for the Chapter 4 algorithms. They encrypt vectors with OpenFHE
CKKS, run the polynomial matrix GROUP BY, lookup JOIN, or RNS digit-mask
construction, decrypt the final slots, and report max absolute error against
an in-process plaintext baseline.
For the RNS experiment, the decomposition step is represented as the data
encoding layout: the client-side encoder materializes encrypted digit columns,
and the online CKKS path evaluates the per-digit Lagrange masks and multiplies
them to recover the target-value mask.

Examples:

```bash
./build/paper_experiments/exp_chap4_group_by --rows 8 --groups 4
./build/paper_experiments/exp_chap4_join --rows 8 --domain 4
./build/paper_experiments/exp_chap4_rns --rows 8 --base 4 --digits 3 --target 15
```

## Chapter 5 TPC-H Query Experiments

The TPC-H targets are compact synthetic reproductions of the query structure,
not a full official TPC-H data-loader benchmark. They are intended for code
review of the encrypted operator pipeline:

- `exp_chap5_query6`: predicate mask times revenue, then rotate-and-sum.
- `exp_chap5_query14`: date predicate mask plus part-type GROUP BY.
- `exp_chap5_query3`: customer/order/lineitem predicate masks, lookup joins,
  and order-key plus ship-priority GROUP BY.
- `exp_chap5_query5`: nation/customer/supplier/order/lineitem lookup chain,
  equality mask, and nation GROUP BY.

Each query accepts:

```bash
--rows N
--seed S
--where-gapmsb
--tfhe-plain-mask
--ckks-mask
--ckks-depth D
--kappa K
--date-bits B
--quantity-bits B
--discount-bits B
--key-bits B
```

The default mode is `--where-gapmsb`. It follows the same high-level shape as
HE3DB's TPC-H Q6 implementation:

1. Encrypt the WHERE operands as TFHE TLWEs.
2. Evaluate each predicate with the Chapter 3 `HomCompare` wrapper, which calls
   the recursive HomMSB/GapMSB path.
3. Convert the TFHE sign result to a 29-bit arithmetic mask and repack the LWE
   masks to OpenFHE CKKS slots through `EvalFHEWtoCKKS`.
4. Run the remaining CKKS filtering, summation, lookup JOIN, equality mask, or
   matrix GROUP BY operator.

`--tfhe-plain-mask` keeps the same TFHE-to-CKKS repacking stage but encrypts
the already-known plaintext predicate mask as TFHE TLWEs. It is useful for
isolating the repacking and CKKS stages. `--ckks-mask` encrypts plaintext masks
directly as CKKS vectors and is only a fast fallback for JOIN/GROUP BY smoke
checks.

TPC-H experiment variables are exposed as command-line parameters instead of
being fixed in the entry points. The query-specific knobs are:

```bash
# Q6
--q6-shipdate-data-min V --q6-shipdate-data-max V
--q6-shipdate-lower V --q6-shipdate-upper V
--q6-discount-lower V --q6-discount-upper V
--q6-quantity-upper V

# Q14
--q14-shipdate-lower V --q14-shipdate-upper V
--q14-type-domain N

# Q3
--q3-key-domain N --q3-priority-domain N --q3-segment V
--q3-orderdate-upper V --q3-shipdate-lower V

# Q5
--q5-key-domain N --q5-nation-domain N --q5-region V
--q5-orderdate-lower V --q5-orderdate-upper V
```

Example smoke commands:

```bash
./build/paper_experiments/exp_chap5_query6 --rows 8 --seed 42
./build/paper_experiments/exp_chap5_query14 --rows 8 --seed 42
./build/paper_experiments/exp_chap5_query3 --rows 8 --seed 42
./build/paper_experiments/exp_chap5_query5 --rows 8 --seed 42
./build/paper_experiments/exp_chap5_query6 --rows 1 --seed 42 --where-gapmsb
./build/paper_experiments/exp_chap5_query6 --rows 4 --seed 42 --tfhe-plain-mask
```

## Security Parameters

`paper_experiments` uses TFHEpp's default 128-bit parameter headers and OpenFHE
`HEStd_128_classic` / `STD128`. The review compile check asserts
`PaperReview::kRequiredSecurityBits == 128`.
