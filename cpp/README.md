# High-performance Fyffe backend

`fyffe_pareto_dp.cpp` is the publication backend. It enumerates capped
capacity PMFs up to a verified feasible cost and removes only first-order
stochastically dominated PMFs at the same layer and exact cost. The pruning is
lossless under the common non-negative future convolutions. The Python kernel
then independently re-evaluates every returned chromosome.

`fyffe_solver.cpp` is the earlier branch-and-bound prototype. It is retained
for audit history but is not the source of the reported 30-cell results.

Example build from this directory:

```text
g++ -O3 -std=c++17 fyffe_pareto_dp.cpp -o fyffe_pareto_dp.exe
fyffe_pareto_dp.exe ..\results\generated\fyffe_pareto_dp.csv
python ..\experiments\verify_fyffe_dp.py
```
