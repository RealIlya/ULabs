# OOP Lab 5 (variant 3)

Implements aggregation-by-reference between `Empiric` (data) and `Estimate`
(M-estimate of shift) using the Observer pattern.

Variant 3 uses the Cauchy weight function:

- `w(z) = c^2 / (c^2 + z^2)`

## Build & run

From repo root:

```bash
g++ -std=c++17 -O2 -I./OOP -o OOP/lab5/demo \
  OOP/lab5/demo.cpp OOP/lab5/Empiric.cpp OOP/lab5/Estimate.cpp OOP/lab5/MainDist.cpp
./OOP/lab5/demo
```
