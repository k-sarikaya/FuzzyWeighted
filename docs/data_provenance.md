# Data provenance and claim boundary

The executable catalogs in `src/fuzzy_wkon/benchmarks.py` are byte-for-value
transcriptions of the V1 implementation. The row-level inventory is in
`benchmark_provenance_tables.md`. Source-paper archiving has begun, but the
implemented catalogs still must be treated as adapted synthetic instances unless
field-by-field reconciliation establishes otherwise.

## Source chain

| Artifact | Immediate source | Upstream publication named by V1 | Verification status |
|---|---|---|---|
| `pharma_adapted` catalog | Discrete instance adapted from Garg and Sharma (2013) | Garg and Sharma (2013), *Reliability-redundancy allocation problem of pharmaceutical plant* | Original continuous RRAP PDF archived under `data/original/`; discrete options adapted for weighted k-out-of-n:G evaluation as detailed in Section 4.1. |
| `fyffe_adapted` catalog | Benchmark instance adapted from Fyffe, Hines and Lee (1968) | Fyffe, Hines and Lee (1968), DOI 10.1109/TR.1968.5217517 | Original Table I archived verbatim in `data/original/fyffe_1968_table_i.csv`. 31 of 42 options match Table I exactly; 11 options adapted as audited in `fyffe_provenance_v6.json` and Section 4.1. |
| `power_rts_adapted` catalog | Generation fleet case study parameterised from the IEEE RTS-79 generating-unit reliability data | IEEE Reliability Test System Task Force (1979), IEEE Trans. PAS-98(6), pp. 2047-2054 | Generating unit classes, capacities (MW), MTTF, MTTR, and FOR values transcribed directly in `data/original/ieee_rts79_units.csv`. Availability probabilities $p_i = \text{MTTF}/(\text{MTTF}+\text{MTTR})$ match the source unit statistics ($1-\text{FOR}$). Peak demand threshold $\tilde{k} = (162, 180, 198)$~MW is a study modeling assumption (regional bay adequacy, not RTS total system peak). Capital expenditure ($M\$) is an illustrative CapEx assumption based on representative $/kW overnight costs. |
| Null option `(0,0,0)` | Already present in every V1 row | No source location demonstrated | Treat as a modeling adaptation (permits leaving a bay or subsystem uncommitted). |
| Capacity TFNs | Symmetric/asymmetric proportional spread rule | Thermodynamic ambient derating | Treat as an engineering modeling assumption justified by seasonal temperature variations. |
| Unit-specific spreads (`derate_low`, `derate_high` in `ieee_rts79_units.csv`) | Hydro 0.10/0; combustion turbine 0.15/0.05; steam units 0.05/0.02 | Hydro: RTS-79 Table 6 (available hydro capacity 100% in quarters 1-2, 90% in quarters 3-4). Thermal: none, RTS-79 Table 4 covers full outages only and leaves partial outages to the user | Hydro taken from the source; thermal values are an assumption of this study (combustion turbines most sensitive to ambient temperature, steam units mildly), with sensitivity sweeps in `rts_unit_specific_v7.json`. Capital costs of the RTS-79 catalog are an illustrative assumption (0.85-1.17 M$/MW); RTS-79 Table 5 gives operating costs only. |
| Fixed reliability targets | Case study constraints | User-specified adequacy target | Treat as a study-specific adequacy requirement ($R_{\text{target}} = 0.95$ for RTS; $0.90$ for Fyffe). |

## Permitted wording

Use "adapted synthetic Pharma/Fyffe-inspired catalog." Do not use "authentic,"
"unchanged," "industrial benchmark reproduction," or "published benchmark data"
unless archived copies of the original source tables are added under
`data/original/` and a field-by-field reconciliation is completed.

For a future faithful claim, record for every numeric field: publication,
table/page/equation, original unit, original value, repository value,
transformation, and verification initials/date. Any unmatched field must remain
explicitly marked `unverified`.
