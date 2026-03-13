# Full-Scale Fission Reactor Benchmark Candidates

Researched 2026-03-13. Selection criteria: full-core models (not pin cells or assemblies),
published geometry specs, reference k-eff/power distributions, tractable on 8GB/10 cores.

| # | Benchmark | Type | Assemblies | Lattice | Reference k-eff | Status |
|---|-----------|------|-----------|---------|----------------|--------|
| 1 | BEAVRS (MIT) | 4-loop W PWR | 193 (17x17) | Square | Measured critical | **FINALIST** |
| 2 | VERA Problem 5 (Watts Bar 1) | W PWR | 193 (17x17) | Square | MC ref ~8 pcm | DISCARDED (too similar to BEAVRS) |
| 3 | NuScale SMR (McSAFER) | SMR PWR | 37 (17x17) | Square | Serpent ref, 44 pcm | **FINALIST** |
| 4 | VVER-1000 MOX (OECD/NEA) | VVER-1000 | 163 hex | Hexagonal | MCU/MCNP ref (6 states) | **FINALIST** |
| 5 | VVER-1000 Temelin | VVER-1000 | 163 hex | Hexagonal | Pin-power dists | DISCARDED (similar to #4) |
| 6 | MC Perf Benchmark | Simplified PWR | 241 (17x17) | Square | ~1.0 | DISCARDED (not realistic) |
| 7 | Peach Bottom 2 BWR | GE BWR/4 | 764 (7x7) | Square | Transient data | DISCARDED (transient-focused) |
| 8 | OPR-1000 | CE PWR | 177 (16x16) | Square | Published | DISCARDED (limited public specs) |
| 9 | PBMR-400 | Pebble-bed HTGR | ~500k pebbles | Random | k=1.06-1.28 | DISCARDED (memory prohibitive) |
| 10 | BN-600 Hybrid Core | SFR | 369 hex | Hexagonal | IAEA-TECDOC-1623 | **FINALIST** |
| 11 | EBR-II Run 138B | SFR | 637 hex | Hexagonal | IRPhE measured | DISCARDED (HEU, restricted data) |
