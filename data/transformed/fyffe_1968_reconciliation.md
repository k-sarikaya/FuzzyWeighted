# Reconciliation: Fyffe Table I vs. `fyffe_adapted` Catalog

## Inputs compared

| Item | Evidence |
|---|---|
| Original values | `data/original/fyffe_hines_lee_1968.pdf`, Table I (p. 68); normalized transcription: `data/original/fyffe_1968_table_i.csv` |
| Adapted catalog values | `src/fuzzy_wkon/benchmarks.py`, `BENCHMARKS["fyffe_adapted"]` |
| Provenance script | `experiments/verify_fyffe_provenance_v6.py` |

## Reconciliation result

The two data structures are structurally adapted. The original problem has 48 active
alternatives (3 or 4 per unit); the adapted catalog has 42 active options (three per unit) plus
14 null options. Tuple equality at a shared unit/alternative index is not a
valid transformation rule: the adapted benchmark alters values and, in several units, replaces
or omits source alternatives. Examples include:

| Unit | Original Table I | Adapted catalog active options | Consequence |
|---|---|---|---|
| 1 | `(0.90,1,3)`, `(0.93,1,4)`, `(0.91,2,2)`, `(0.95,2,5)` | `(0.90,1,3)`, `(0.93,1,4)`, `(0.94,2,5)` | Drops two source options and introduces `(0.94,2,5)`. |
| 2 | `(0.95,2,8)`, `(0.94,1,10)`, `(0.93,1,9)` | `(0.95,2,8)`, `(0.94,1,10)`, `(0.96,3,11)` | First two values match; third is changed. |
| 6 | `(0.99,3,5)`, `(0.98,3,4)`, `(0.97,2,5)`, `(0.96,2,4)` | `(0.99,3,5)`, `(0.98,3,4)`, `(0.99,4,6)` | Retains two source values, omits two, and adds one changed tuple. |
| 14 | `(0.90,4,6)`, `(0.92,4,7)`, `(0.95,5,6)`, `(0.99,6,9)` | `(0.90,4,6)`, `(0.92,4,7)`, `(0.95,5,8)` | Retains two source values, alters the third, and omits the fourth. |

More fundamentally, the source tuple is `(unit reliability, unit cost, unit
weight)` and allows an integer redundancy count. The adapted benchmark treats its third tuple
field as a capacity, includes a null option, applies capacity fuzzification,
and solves a discrete-choice weighted $k$-out-of-$n$:G problem subject to a
prescribed reliability target. These changes are intentional model adaptations
documented in Section~4.1 of the manuscript.

## Approved claim language

The original Table I is archived and traceable. The benchmark catalog represents an
**adapted synthetic Fyffe-inspired benchmark**; as disclosed in the manuscript, it adapts
the physical topology and catalog options to the weighted $k$-out-of-$n$:G formulation.
Source recovery enables this provenance to be transparently documented and independently verified.
