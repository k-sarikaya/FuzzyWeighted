# Fyffe, Hines and Lee (1968) source manifest

## Archived primary source

- File: `fyffe_hines_lee_1968.pdf`
- SHA256: `A2AC7A8ED5A632E77158E5782645C8F9ED1AA2DF32BAC5F5A33DEB6788CCFA28`
- Citation: D. E. Fyffe, W. W. Hines and N. K. Lee, "System Reliability
  Allocation and a Computational Algorithm", *IEEE Transactions on
  Reliability*, R-17(2), 64--69, 1968, DOI
  [10.1109/TR.1968.5217517](https://doi.org/10.1109/TR.1968.5217517).
- Source table: Table I, p. 68, "Numerical Input Data for a Reliability
  Allocation Problem".
- Archived: 2026-08-12.

The NASA NTRS record `19680054179` remains useful as bibliographic metadata,
but it does not itself expose a downloadable file. The article PDF above was
subsequently obtained from an IEEE-hosted copy and is now archived in this
repository.

## What the source actually specifies

The paper's numerical example is a **14-unit series reliability-allocation
problem**. For each functional unit, Table I gives 3--4 design alternatives;
each alternative has unit reliability `R`, unit cost `C`, and unit weight `W`.
The source model selects a design alternative and an integer number of
identical parallel units per functional unit, maximizes series-system
reliability, and imposes total cost and weight constraints. Its illustrated
case uses cost allowance 130 and weight allowance 170 (p. 69).

`fyffe_1968_table_i.csv` is a direct transcription of Table I. It is a
source-data artifact, not an input to the V2 solver.

## Claim boundary after source recovery

Recovery of Table I removes the former "source table unavailable" limitation,
but it does **not** make the V2 executable `fyffe_adapted` catalog a faithful
Fyffe benchmark. The V2 catalog adds a null option, retains only three active
options per unit, changes several numerical tuples, reinterprets source weight
as operational capacity, replaces the source redundancy decision with a
finite-choice weighted-system problem, and adds TFNs and a fixed reliability
target. The field-by-field reconciliation is recorded in
`../transformed/fyffe_1968_reconciliation.md`.

Permitted wording for V2 is therefore "Fyffe-inspired adapted synthetic
catalog, reconciled against Table I of Fyffe et al. (1968)". A faithful
reproduction claim remains prohibited for the current executable catalog.
