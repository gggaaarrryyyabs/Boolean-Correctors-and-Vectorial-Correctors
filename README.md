The Boolean and vectorial Boolean corrector examples of [Constructions of Boolean Correctors with Optimal Nonlinearity and Vectorial Correctors Beyond Optimal Linear Correction Order](https://eprint.iacr.org/2026/2812) in Cryptology ePrint Archive.


## Contents

`/data` in this repo contains $14$ Boolean correctors of correction order $n−2$ for $n = 3,\dots,16$, four Boolean correctors of correction order $n−3$ for $n = 5,...,8$, and five vectorial correctors. 
``verify_structure_numpy.py`` and ``verify_structure.py`` are verifying Python files, where the first file requires NumPy.
The verification results will be `/verification/verification_independent.json` and `/verification/verification_numpy.json`, respectively.

| Content in /data | File or file prefix |
| --- | --- |
| $(n-2)$-correctors ($3\le n \le 10$) | `best_corrector_n*.txt` |
| $(n-2)$-correctors ($n=11,12,13$) | `rotation_quota_C10_n11.txt`, `rotation_quota_C11_n12.txt`, `rotation_quota_C12_n13.txt` |
| $(n-2)$-correctors ($n=14,15,16$) | `layer_quota_local_n14.txt`, `layer_quota_local_n15.txt`, `layer_quota_local_n16.txt` |
| $(n-3)$-correctors ($n=5,\dots,8$) | `corrector_n5_c2_bound_truth_table.txt`, `corrector_n6_c3_bound_truth_table.txt`, `corrector_n7_c4_candidate_truth_table.txt`, `corrector_n8_c5_M24_truth_table.txt` |
| Gold vectorial correctors | `gold_vector_seed_n5_m2_c1_outputs.txt`, `gold_vector_seed_n7_m2_c2_outputs.txt` |
| Welch $(5,3,1)$-correctors | `welch_vector_seed_n5_m3_c1_outputs.txt` |
| $(9,3,4)$-correctors | `vector_corrector_n9_m3_c4_L40_M48_outputs.txt` |
| $(18,2,13)$-correctors | `vector_corrector_n18_m2_c13_M1728_outputs.bin` |

## Encoding

Inputs occur in increasing integer order $x=\sum_i x_i 2^i$, thus $x_0$ is the least-significant input bit. 
A Boolean text file is a $0/1$ truth-table string, wrapped at $64$ characters. 
A vectorial text file has one decimal output per line, with $y = \sum_i y_i 2^i$. 

The $n=18$ binary file has exactly $262,144$ bytes, one unsigned byte per input. 
Every byte is $0, 1, 2$, or $3$. 

## Verified parameters

All listed resiliency values and correction orders are exact. 
Vectorial nonlinearity and resiliency are minimum over all nonzero components, maximum Walsh magnitude is the maximum over all components, and algebraic degree is the maximum algebraic degree over all coordinate functions.

| File | n | m | c | M | NL | res | deg |
|---|---:|---:|---:|---:|---:|---:|---:|
| `best_corrector_n3.txt`  | 3  | 1 | 1 | 4   | 2   | 0 | 2 |
| `best_corrector_n4.txt`  | 4  | 1 | 2 | 8   | 4   | 0 | 2 |
| `best_corrector_n5.txt`  | 5  | 1 | 3 | 8   | 12  | 0 | 3 |
| `best_corrector_n6.txt`  | 6  | 1 | 4 | 12  | 26  | 0 | 5 |
| `best_corrector_n7.txt`  | 7  | 1 | 5 | 20  | 54  | 0 | 6 |
| `best_corrector_n8.txt`  | 8  | 1 | 6 | 32  | 112 | 0 | 6 |
| `best_corrector_n9.txt`  | 9  | 1 | 7 | 60  | 226 | 0 | 8 |
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

## Verification

From the extracted package directory:

```bash
python verify_structure_numpy.py
python verify_structure.py
```

The first program requires NumPy, the second uses only the Python standard library. 
Both exit with status 0 only when all 23 files pass. 
Both recompute full Walsh spectra and algebraic degrees directly from the truth tables. 


## Environment

The two Python files require the following environment:
```
"python": "3.13.5"
"numpy": "2.3.5"
```
