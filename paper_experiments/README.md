# Paper Experiment Runbook

This directory is for reviewable paper-facing experiments. The executables do
not print fixed manuscript tables. Each target constructs input data, encrypts
the columns or masks, runs the homomorphic operator under test, decrypts only
the final result, and checks it against a plaintext baseline computed in the
same executable.

## Build

```bash
cmake -S . -B build -DENABLE_TEST=ON
cmake --build build --target paper_algorithms paper_native_baselines \
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
2. Plaintext case generation for unsigned operands. A `bits` value means the
   input operands are in `[0, 2^bits-1]`; the encrypted signed difference uses
   one extra sign bit, matching Chapter 3's random-input experiment.
3. TFHE lvl2 TLWE encryption of `lhs` and `rhs`.
4. `PaperReview::Chapter3HomCompare`, which dispatches to the `metapbs2`
   recursive HomMSB/GapMSB comparison path.
5. lvl1 sign decryption and row-by-row comparison with the plaintext predicate.
   Relation, equality, and inequality predicates are selected with
   `--predicate relation|equality|inequality`; equality and inequality compute
   both `lhs > rhs` and `lhs < rhs` under TFHE and combine the encrypted bits
   before the final decryption.

`exp_chap3_k_impact` follows the Chapter 3 zeroed-bit experiment directly: it
encrypts boundary-neighborhood messages, runs `GapMSB_k` for each configured
zeroed bit `k`, and reports `w_k`, the first-round pruning period, pruning
statistics, latency, and the decrypted MSB check.

Examples:

```bash
./build/paper_experiments/exp_chap3_accuracy --bits-list 4,8,16,32 --trials 100
./build/paper_experiments/exp_chap3_compare_test --bits 16 --predicate equality --trials 20 --seed 7
./build/paper_experiments/exp_chap3_boundary --bits 16 --threshold 32768 --radius 4 --kappa 5
./build/paper_experiments/exp_chap3_k_impact --k-list 10,9,8,7,1 --radius 4 --trials 100
./build/paper_experiments/exp_chap3_compare_test --bits 16 --predicate inequality --output results/ch3_neq.csv
```

## Chapter 4 CKKS Operators

`exp_chap4_group_by`, `exp_chap4_join`, and `exp_chap4_rns` are operator
experiments for the Chapter 4 algorithms. They encrypt vectors with OpenFHE
CKKS, run the polynomial matrix GROUP BY, lookup JOIN, or RNS digit-mask
construction, decrypt the final slots, and report correctness/accuracy plus
absolute-error statistics against an in-process plaintext baseline. The
decrypted-slot correctness threshold is fixed in the code at `1e-6`; it is not
an experiment parameter. Each executable prints the same variables used in the
Chapter 4 tables: problem size, domain or group count, RNS degree/configuration,
correctness/accuracy, and absolute-error metrics.
For the RNS experiment, the decomposition step is represented as the data
encoding layout: the client-side encoder materializes encrypted digit columns,
and the online CKKS path evaluates the per-digit Lagrange masks and multiplies
them to recover the target-value mask.

Examples:

```bash
./build/paper_experiments/exp_chap4_group_by --rows 4096 --group-bits 4 --system ours
./build/paper_experiments/exp_chap4_group_by --rows 4096 --group-bits 4 --system he3db
./build/paper_experiments/exp_chap4_group_by --rows 4096 --group-bits 4 --system arcedb
./build/paper_experiments/exp_chap4_join --rows 4096 --key-domain 64 --seed 42 --system ours
./build/paper_experiments/exp_chap4_join --rows 4096 --key-domain 64 --seed 42 --system he3db
./build/paper_experiments/exp_chap4_join --rows 4096 --key-domain 64 --seed 42 --system arcedb
./build/paper_experiments/exp_chap4_rns --rows 4096 --domain-size 256 --base 4 --target 15 --output results/ch4_rns.csv
```

`--system he3db` and `--system arcedb` run enumerated equality baselines using
vendored native comparison operators under
`paper_experiments/baselines/he3db_native` and
`paper_experiments/baselines/arcedb_native`. They do not claim that HE3DB or
ArcEDB has native GROUP BY or JOIN. The GROUP BY/JOIN logic is the paper's
enumerated equality baseline: run the native equality comparison once per
candidate group or join key, then aggregate or backfill the matching payload.
The output reports `enumerated_groups`, `enumerated_keys`, `equality_calls`,
`time_source=measured_native_operator`, and the same accuracy/error fields as
the CKKS operator path.

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
--row-exp E
--seed S
--output PATH
--ckks-depth D
--kappa K
--date-bits B
--quantity-bits B
--discount-bits B
--key-bits B
--system ours|he3db|arcedb|all
```

Chapter 5 uses a single online mask path. It follows the same high-level shape
as HE3DB's TPC-H Q6 implementation:

1. Encrypt the WHERE operands as TFHE TLWEs.
2. Evaluate each predicate with the Chapter 3 `HomCompare` wrapper, which calls
   the recursive HomMSB/GapMSB path.
3. Convert the TFHE sign result to a 29-bit arithmetic mask and repack the LWE
   masks to OpenFHE CKKS slots through `EvalFHEWtoCKKS`.
4. Run the remaining CKKS filtering, summation, lookup JOIN, equality mask, or
   matrix GROUP BY operator.

The Chapter 5 timing summary follows HE3DB's `test/tpch_q6.cpp` split:
`filter_time_ms` contains only the online TFHE predicate comparison loop.
`aggregation_time_ms` includes TFHE-to-CKKS repacking plus CKKS mask products,
SUM rotation, JOIN, equality-mask, and GROUP BY work. Key generation, synthetic
data generation, CKKS column encryption, and setup stages are printed separately
under `stage_category=setup` and are not included in
`total_query_time_ms = filter_time_ms + aggregation_time_ms`.

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
--q14-shipdate-data-min V --q14-shipdate-data-max V
--q14-type-domain N --q14-promo-type V

# Q3
--q3-key-domain N --q3-priority-domain N --q3-segment V
--q3-segment-domain N
--q3-orderdate-base V --q3-orderdate-step V
--q3-shipdate-data-min V --q3-shipdate-data-max V
--q3-orderdate-upper V --q3-shipdate-lower V

# Q5
--q5-key-domain N --q5-nation-domain N --q5-region-domain N --q5-region V
--q5-orderdate-base V --q5-orderdate-step V
--q5-orderdate-lower V --q5-orderdate-upper V
```

Example smoke commands:

```bash
./build/paper_experiments/exp_chap5_query6 --row-exp 10 --seed 42
./build/paper_experiments/exp_chap5_query14 --row-exp 10 --seed 42 --q14-type-domain 4
./build/paper_experiments/exp_chap5_query3 --row-exp 10 --seed 42 --q3-key-domain 16 --q3-priority-domain 4 --system ours
./build/paper_experiments/exp_chap5_query3 --row-exp 10 --seed 42 --q3-key-domain 16 --q3-priority-domain 4 --system he3db
./build/paper_experiments/exp_chap5_query3 --row-exp 10 --seed 42 --q3-key-domain 16 --q3-priority-domain 4 --system arcedb
./build/paper_experiments/exp_chap5_query5 --row-exp 10 --seed 42 --q5-key-domain 16 --q5-nation-domain 8 --system ours
./build/paper_experiments/exp_chap5_query5 --row-exp 10 --seed 42 --q5-key-domain 16 --q5-nation-domain 8 --system he3db
./build/paper_experiments/exp_chap5_query5 --row-exp 10 --seed 42 --q5-key-domain 16 --q5-nation-domain 8 --system arcedb
./build/paper_experiments/exp_chap5_query6 --rows 1 --seed 42 --output results/q6.csv
```

For Q3 and Q5, `--system he3db` and `--system arcedb` execute the Chapter 5
functional-equivalent enumerated JOIN baseline. WHERE/JOIN/group masks are
generated by the vendored native comparison operators, while the query plan
enumerates candidate keys because HE3DB and ArcEDB do not provide native JOIN
or GROUP BY operators. `--system all` prints the HE3DB baseline, ArcEDB
baseline, and then the current CKKS `ours` path.

## Security Parameters

`paper_experiments` uses TFHEpp's default 128-bit parameter headers and OpenFHE
`HEStd_128_classic` / `STD128`. The review compile check asserts
`PaperReview::kRequiredSecurityBits == 128`.
