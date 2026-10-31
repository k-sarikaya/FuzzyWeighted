# Source-level benchmark provenance tables

Each option is written `(reliability, cost, nominal capacity)`. `O0` is the
null option. The immediate source is the V1 code line shown. The upstream
Pharma paper has now been archived, and the reconciliation file
`data/transformed/pharma_original_reconciliation.csv` shows that the V2 catalog
is structurally different from the original continuous RRAP. The upstream Fyffe
Table I is now archived and transcribed under `data/original/`, and its
comparison with V2 is recorded in `data/transformed/fyffe_1968_reconciliation.md`.
Consequently every V2 row remains an **adaptation**, not a faithful reproduction.

## Pharma-adapted catalog

Global fields: `N=10` (V1 line 50), `K_tfn=(34,40,46)` (line 51),
`R_target=0.85` (line 55). The TFN and fixed target have no demonstrated
location in Garg and Sharma (2013). The archived source article instead defines
a continuous reliability-redundancy allocation problem with component
reliabilities `r_i`, redundancy counts `n_i`, and nonlinear cost, weight and
volume constraints.

| Subsystem | O0 | O1 | O2 | Immediate source | Status |
|---|---|---|---|---|---|
| S1 | (0,0,0) | (0.1234,10,9) | (0.1500,12,10) | V1 line 57 | Unverified adaptation |
| S2 | (0,0,0) | (0.8838,20,7) | (0.8900,22,8) | V1 line 58 | Unverified adaptation |
| S3 | (0,0,0) | (0.8627,15,5) | (0.8800,18,6) | V1 line 59 | Unverified adaptation |
| S4 | (0,0,0) | (0.8666,12,9) | (0.8900,15,10) | V1 line 60 | Unverified adaptation |
| S5 | (0,0,0) | (0.4509,8,9) | (0.5000,10,10) | V1 line 61 | Unverified adaptation |
| S6 | (0,0,0) | (0.8223,25,10) | (0.8500,30,12) | V1 line 62 | Unverified adaptation |
| S7 | (0,0,0) | (0.7105,18,6) | (0.7500,20,8) | V1 line 63 | Unverified adaptation |
| S8 | (0,0,0) | (0.7930,14,5) | (0.8100,16,6) | V1 line 64 | Unverified adaptation |
| S9 | (0,0,0) | (0.6990,12,8) | (0.7200,15,9) | V1 line 65 | Unverified adaptation |
| S10 | (0,0,0) | (0.8593,16,6) | (0.8800,18,7) | V1 line 66 | Unverified adaptation |

## Fyffe-adapted catalog

Global fields: `N=14` (V1 line 70), `K_tfn=(40,50,60)` (line 71),
`R_target=0.90` (line 75). The TFN, fixed target, null options, and the use of
the third tuple field as operational capacity differ from Fyffe, Hines and Lee
(1968). The original Table I is now archived; see
`data/transformed/fyffe_1968_reconciliation.md` for field-level examples.

| Subsystem | O0 | O1 | O2 | O3 | Immediate source | Status |
|---|---|---|---|---|---|---|
| S1 | (0,0,0) | (0.90,1,3) | (0.93,1,4) | (0.94,2,5) | V1 line 77 | Unverified adaptation |
| S2 | (0,0,0) | (0.95,2,8) | (0.94,1,10) | (0.96,3,11) | V1 line 78 | Unverified adaptation |
| S3 | (0,0,0) | (0.85,2,7) | (0.90,3,5) | (0.92,4,6) | V1 line 79 | Unverified adaptation |
| S4 | (0,0,0) | (0.83,3,5) | (0.87,4,6) | (0.89,5,7) | V1 line 80 | Unverified adaptation |
| S5 | (0,0,0) | (0.94,2,4) | (0.93,2,3) | (0.95,3,5) | V1 line 81 | Unverified adaptation |
| S6 | (0,0,0) | (0.99,3,5) | (0.98,3,4) | (0.99,4,6) | V1 line 82 | Unverified adaptation |
| S7 | (0,0,0) | (0.91,4,7) | (0.92,4,8) | (0.94,5,9) | V1 line 83 | Unverified adaptation |
| S8 | (0,0,0) | (0.81,3,4) | (0.90,5,7) | (0.92,6,8) | V1 line 84 | Unverified adaptation |
| S9 | (0,0,0) | (0.97,2,8) | (0.99,3,9) | (0.99,4,10) | V1 line 85 | Unverified adaptation |
| S10 | (0,0,0) | (0.83,4,6) | (0.85,4,5) | (0.88,5,6) | V1 line 86 | Unverified adaptation |
| S11 | (0,0,0) | (0.94,3,5) | (0.95,4,6) | (0.96,5,7) | V1 line 87 | Unverified adaptation |
| S12 | (0,0,0) | (0.79,2,4) | (0.82,3,5) | (0.85,4,6) | V1 line 88 | Unverified adaptation |
| S13 | (0,0,0) | (0.98,2,5) | (0.99,3,5) | (0.99,4,6) | V1 line 89 | Unverified adaptation |
| S14 | (0,0,0) | (0.90,4,6) | (0.92,4,7) | (0.95,5,8) | V1 line 90 | Unverified adaptation |

## Transformation register

| Transformation | Rule | Provenance status |
|---|---|---|
| Capacity fuzzification | `(0.9w, w, 1.1w)` in the default scenario | V2 modeling assumption; not attributed to either source paper |
| Capacity scalarization | lower alpha-cut endpoint plus `gamma` times interval width | Declared Hurwicz rule |
| Threshold scalarization | upper alpha-cut endpoint minus `gamma` times interval width | Declared conservative-opposite Hurwicz rule |
| Integerization | exact common-denominator rational scaling | Computational transformation, value-preserving |
| Objective | minimum feasible cost, then maximum reliability at that cost | V2 study objective; not asserted as the source papers' objective |
