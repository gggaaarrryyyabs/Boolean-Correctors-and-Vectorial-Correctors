# Verified appendix corrector certificates

**All 23 truth tables pass both independent verification programs.**

These are newly generated witnesses matching the parameters in the supplied `one-column.tex`, not recovered copies of the original supplementary data. The uploaded source labels the filename table `tab:certificate_files`; the request referred to `tab:example_files`. This package follows the active filenames in that table, excluding commented-out entries.

## Contents

`data/` contains 14 Boolean witnesses of exact order n−2 (n = 3,…,16), four Boolean witnesses of exact order n−3 (n = 5,…,8), and five vectorial witnesses. The two verification scripts use the filenames listed in the appendix. `manifest.json` records the expected parameters and file formats. `SHA256SUMS` records the hashes of all truth tables. Detailed component spectra, layer sums, and preimage-layer counts are in `verification/verification_numpy.json`; an independently implemented check is in `verification/verification_independent.json`.

## Encoding: input and output order

Inputs occur in increasing integer order x = sum_i x_i 2^i. Thus x_0 is the least-significant input bit. A Boolean text file is an ASCII 0/1 truth-table string, wrapped at 64 characters; whitespace is not part of the table. A vectorial text file has one decimal output label per line, with y = sum_i y_i 2^i. Component lambda is parity(lambda & y).

The n=18 binary file has exactly 262,144 bytes, one unsigned byte per input. Every byte is 0, 1, 2, or 3. It is NOT a packed two-bit file, a hexadecimal string, or a text file.

## Verified parameters

All listed resiliency values are exact; all listed correction orders are exact. Vectorial nonlinearity and resiliency are minima over all nonzero components, maximum Walsh magnitude is the maximum over all components and frequencies, and algebraic degree is the maximum coordinate/component degree.

| File | n | m | c | M | NL | res | degree |
|---|---:|---:|---:|---:|---:|---:|---:|
| `best_corrector_n3.txt` | 3 | 1 | 1 | 4 | 2 | 0 | 2 |
| `best_corrector_n4.txt` | 4 | 1 | 2 | 8 | 4 | 0 | 2 |
| `best_corrector_n5.txt` | 5 | 1 | 3 | 8 | 12 | 0 | 3 |
| `best_corrector_n6.txt` | 6 | 1 | 4 | 12 | 26 | 0 | 5 |
| `best_corrector_n7.txt` | 7 | 1 | 5 | 20 | 54 | 0 | 6 |
| `best_corrector_n8.txt` | 8 | 1 | 6 | 32 | 112 | 0 | 6 |
| `best_corrector_n9.txt` | 9 | 1 | 7 | 60 | 226 | 0 | 8 |
| `best_corrector_n10.txt` | 10 | 1 | 8 | 104 | 460 | 0 | 9 |
| `rotation_quota_C10_n11.txt` | 11 | 1 | 9 | 188 | 930 | 0 | 10 |
| `rotation_quota_C11_n12.txt` | 12 | 1 | 10 | 344 | 1876 | 0 | 10 |
| `rotation_quota_C12_n13.txt` | 13 | 1 | 11 | 632 | 3780 | 0 | 11 |
| `layer_quota_local_n14.txt` | 14 | 1 | 12 | 1172 | 7606 | 0 | 13 |
| `layer_quota_local_n15.txt` | 15 | 1 | 13 | 2188 | 15290 | 0 | 14 |
| `layer_quota_local_n16.txt` | 16 | 1 | 14 | 4096 | 30720 | 0 | 14 |
| `corrector_n5_c2_bound_truth_table.txt` | 5 | 1 | 2 | 8 | 12 | 0 | 3 |
| `corrector_n6_c3_bound_truth_table.txt` | 6 | 1 | 3 | 12 | 26 | 0 | 5 |
| `corrector_n7_c4_candidate_truth_table.txt` | 7 | 1 | 4 | 16 | 56 | 0 | 4 |
| `corrector_n8_c5_M24_truth_table.txt` | 8 | 1 | 5 | 24 | 116 | 0 | 6 |
| `gold_vector_seed_n5_m2_c1_outputs.txt` | 5 | 2 | 1 | 8 | 12 | 0 | 2 |
| `gold_vector_seed_n7_m2_c2_outputs.txt` | 7 | 2 | 2 | 16 | 56 | 0 | 2 |
| `welch_vector_seed_n5_m3_c1_outputs.txt` | 5 | 3 | 1 | 8 | 12 | 0 | 3 |
| `vector_corrector_n9_m3_c4_L40_M48_outputs.txt` | 9 | 3 | 4 | 48 | 232 | 0 | 7 |
| `vector_corrector_n18_m2_c13_M1728_outputs.bin` | 18 | 2 | 13 | 1728 | 130208 | 0 | 17 |

## Component-polynomial checks

For each n−2 Boolean witness, the exact Walsh-layer polynomial is 2^n z^(n−1). For each n−3 Boolean witness it is 2^n z^(n−2).

For the (9,3,4) witness the exact component polynomials, using the output-bit convention above, are:

```text
lambda = 1,5,7: 512 (2 z^5 - z^7)
lambda = 2,4,6: 512 z^6
lambda = 3:     512 (-z^5 + 2 z^7)
```

For the (18,2,13) witness, every nonzero component has polynomial 2^18 (2 z^14 - z^16). Its verified aggregate parameters are M=1728, NL=130208, res=0, degree=17. The smaller-degree intermediate found during generation is NOT included in `data/`.

## Verification

From the extracted package directory:

```bash
python verify_revision_structure.py
python verify_balanced_field_orbit.py
sha256sum -c SHA256SUMS
```

The first program requires NumPy; the second uses only the Python standard library. Both exit with status 0 only when all 23 files pass. Both recompute full Walsh spectra and algebraic degrees directly from the truth tables. The second program separately computes Walsh-layer sums from preimage counts using exact Krawtchouk coefficients, and checks the recorded finite-field realizations of the Gold/Welch examples. No external SAT/MILP solver is needed to verify the witnesses.

## Construction and provenance

The n≤10 Boolean witnesses and the n−3 witnesses were found by layer-quota-preserving local search. For n=11,…,16, random layer assignments were followed by exact coordinate-incidence balancing. The Gold/Welch examples use explicitly recorded finite fields, input bases and trace coefficients in `verification/field_constructions.json`. The n=9 witness uses its prescribed joint profile and degree-preserving paired swaps. The n=18 witness uses its two-term profile, local swaps, and a final search maintaining exact degree 17.

The names `rotation_quota_C10_n11.txt`, `rotation_quota_C11_n12.txt`, and `rotation_quota_C12_n13.txt` are retained only to match the appendix. These newly generated witnesses are not asserted to possess rotational invariance. Likewise, the otherwise unexplained `L40` token is retained in the n=9 filename without asserting a separate L40 property. The code and metadata, not historical filenames, describe how these witnesses were produced.

All search programs use fixed seeds recorded in the metadata. Source files are supplied under `src/`; compiled executables and intermediate truth tables are not needed for verification. To regenerate the entire package with Python, NumPy and a C++17 compiler:

```bash
python regenerate.py
```

Regeneration has explicit search budgets and raises an error on failure rather than substituting an unverified witness. Random-library or compiler changes may yield a different witness; every replacement must still pass both verifiers. Search failure is never treated as proof of infeasibility.

## Changes needed in the manuscript

These files support the numerical existence claims and match every parameter requested by the active filename table. Their verification is not a replay of the previously claimed CP-SAT, exhaustive orbit, or Windows-based computations. Replace those method/hardware claims where they purport to describe these particular data. The supplied `appendix_certificate_files.tex` provides replacement data/provenance text and a filename table. Analytical lower bounds remain the basis for any optimality claim: the scripts verify witnesses, not the completeness of a search.

The n=8 and n=16 files here match the search-table parameters res=0, degrees 6 and 14. They are not the separate self-addressing theorem instances, whose resiliency and degree differ. No truth table of the enormous lifted examples (for example n=72 or n=837) is claimed to be included.
