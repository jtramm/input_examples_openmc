# Fission Reactor Benchmark Candidates for OpenMC Implementation

This document catalogs physically realistic fission reactor benchmark problems that could
be implemented as OpenMC Monte Carlo input models. Each entry has been verified to have
publicly accessible specifications with sufficient geometry, material, and validation data.

---

## 1. VERA Core Physics Benchmark -- Problem 2 (2D HZP Fuel Lattice)

- **Type:** PWR fuel assembly lattice (critical configuration)
- **Description:** A 2D hot-zero-power beginning-of-cycle fuel lattice based on the
  Westinghouse 17x17 PWR design for Watts Bar Unit 1. The benchmark defines a single
  fuel assembly with 264 fuel rods, 24 guide tubes, and 1 instrument tube, with
  discrete Pyrex burnable poison rods. Assembly enrichments of 2.11%, 2.619%, and 3.10%
  are specified. Problems 1-3 range from a single pin cell to a 3D assembly.
- **Specification URL:** https://corephysics.com/docs/CASL-U-2012-0131-004.pdf
- **Reference results:** Monte Carlo reference solutions for k-eff, pin power
  distributions; VERA results with RMSE < 0.04% for eigenvalues and < 1.3% for pin
  powers have been demonstrated with OpenMC.
- **Complexity:** ~200-500 surfaces for a single assembly; 17x17 pin lattice with
  fuel/gap/clad/moderator detail. Problem 1 (pin cell) is very simple; Problem 2
  (assembly) is moderate; Problem 3 (3D assembly) adds axial complexity.
- **Notes:** Already validated with OpenMC (see S1738573323002632). The CASL-U-2012-0131-004
  document is freely available. Excellent progression from simple to complex. The
  specification PDF appears to be image-based, but the problem geometry is widely
  reproduced in published papers.

---

## 2. BEAVRS (Benchmark for Evaluation And Validation of Reactor Simulations)

- **Type:** Full-core PWR depletion benchmark (commercial plant)
- **Description:** A two-cycle full-core PWR depletion benchmark based on a 4-loop
  Westinghouse plant with 193 optimized fuel assemblies (OFAs) of 17x17 lattice for
  3411 MWth rated power. Provides detailed geometry of fuel assemblies, burnable
  absorbers, in-core fission detectors, core loading/reloading patterns, and numerous
  in-vessel components including baffle and barrel.
- **Specification URL:** https://github.com/mit-crpg/BEAVRS (includes OpenMC model)
- **Reference results:** Measured HZP physics test data, boron letdown curves, and 3D
  in-core flux maps from 58 instrumented assemblies for two full cycles.
- **Complexity:** Full core is very complex (~50,000+ surfaces); however, a single
  assembly or a small subset of assemblies is tractable (~500-2000 surfaces).
- **Notes:** OpenMC Python API model already exists in the repository. Best used by
  extracting a single assembly or small assembly cluster rather than modeling the full
  core. Extremely well-documented and freely available on GitHub.

---

## 3. NuScale-Like SMR Benchmark (McSAFER Project)

- **Type:** Small modular reactor full core
- **Description:** A NuScale-like PWR core with 37 fuel assemblies of 7 different types,
  160 MWth expected power output. Uses standard 17x17 lattice with 264 fuel rods, 24
  guide tubes per assembly, 1.2598 cm fuel rod pitch, 21.5036 cm assembly pitch, and
  200 cm active fuel length. Enrichments of 1.6%, 2.4%, and 3.1% UO2.
- **Specification URL:** https://www.mdpi.com/2673-4362/6/4/44 (OpenMC implementation)
  and https://www.sciencedirect.com/science/article/pii/S1738573323002930 (benchmark definition)
- **Reference results:** Serpent reference solutions for k-eff, control rod worth,
  radial/axial power distributions. OpenMC vs. Serpent k-eff differences of ~355 pcm.
- **Complexity:** Full 37-assembly core is ~5,000-10,000 surfaces; a single assembly is
  ~300-500 surfaces. Boron-free design adds modeling interest.
- **Notes:** Already implemented in OpenMC v0.15.0. The ECP ExaSMR benchmark
  (https://github.com/mit-crpg/ecp-benchmarks) provides a similar NuScale-like model
  with OpenMC Python API scripts.

---

## 4. VENUS-2 MOX Core Experiment

- **Type:** Critical experiment with mixed UO2/MOX fuel
- **Description:** An OECD/NEA benchmark based on the Belgian VENUS-2 reactor, featuring
  a mixed core of 3.3 wt% UO2, 4.0 wt% UO2, and 2.0/2.7 wt% MOX fuel pins. The core
  includes a baffle, reflector, barrel, water gap, and neutron pad regions. Provides a
  realistic representation of PWR fuel management with mixed oxide fuel.
- **Specification URL:** https://www.oecd-nea.org/upload/docs/application/pdf/2019-12/nsc-doc2005-22.pdf
  (OECD/NEA benchmark report, NEA/NSC/DOC(2005)22)
- **Reference results:** Experimental k-eff, axial fission rate distributions, equivalent
  fission fluxes. Multiple code comparison results (TRIPOLI-4, Serpent2, MCNP6.1).
- **Complexity:** ~1,000-3,000 surfaces for the 2D core model with three fuel types.
- **Notes:** Already validated with OpenMC (see Springer article s13369-023-08684-x).
  The benchmark includes both a core physics experiment and a dosimetry response
  experiment. Specification PDF is freely available from OECD/NEA.

---

## 5. MSRE (Molten Salt Reactor Experiment) Benchmark

- **Type:** Research reactor (molten salt, graphite moderated)
- **Description:** Benchmark of the ORNL Molten Salt Reactor Experiment first criticality,
  cataloged in IRPhEP as MSRE-MSR-EXP-001. The reactor vessel has 147.32 cm inner
  diameter and ~238.76 cm height, containing vertical graphite stringers that create
  fuel salt channels. Zero-power experiments with U-235 molten salt fuel from June 1965.
- **Specification URL:** https://www.frontiersin.org/journals/nuclear-engineering/articles/10.3389/fnuen.2024.1385478/full
  (open access; contains CSG and CAD model details for OpenMC)
- **Reference results:** Experimental k-eff = 0.99978 +/- 0.00420; flux distributions.
  Serpent and OpenMC CSG models agree within 10 pcm; spatial/energy flux distributions
  agree within 0.1%.
- **Complexity:** ~500-2,000 surfaces for the CSG model. The graphite stringer array and
  vessel internals provide moderate geometric complexity.
- **Notes:** Both CSG and CAD models have been implemented in OpenMC. The CAD model
  captures more geometric detail and gives ~1% lower k-eff, suggesting physical
  importance of detailed geometry. Unique non-LWR benchmark.

---

## 6. HTR-10 Pebble Bed Reactor

- **Type:** Research reactor (high-temperature gas-cooled, pebble bed)
- **Description:** The Chinese HTR-10 is a 10 MWth pebble bed HTGR with TRISO-coated
  particle fuel. Each fuel pebble contains 8,335 TRISO UO2 kernel coated particles
  at 17% enrichment dispersed in a graphite matrix. The core contains tens of thousands
  of randomly distributed fuel and moderator pebbles. Helium cooled, graphite moderated.
- **Specification URL:** https://www.researchgate.net/publication/348962701_Criticality_Analysis_of_HTR-10_using_an_Open_Source_Monte_Carlo_code_OpenMC
  (open access paper with OpenMC model details)
- **Reference results:** Critical height measurements (experimental: ~123-126 cm depending
  on configuration); VSOP, MCNP, and Serpent reference k-eff values.
- **Complexity:** TRISO particle double heterogeneity requires OpenMC's random packing
  features; ~500-2,000 surfaces depending on packing model fidelity.
- **Notes:** Demonstrates OpenMC's TRISO/pebble modeling capabilities. The double
  heterogeneity (TRISO in pebbles, pebbles in core) is a distinctive modeling challenge.
  Well-suited as an advanced reactor benchmark.

---

## 7. KRITZ-2 Critical Experiments (LEU and MOX Lattices)

- **Type:** Critical experiment (fuel rod lattice at varying temperatures)
- **Description:** Thermal critical experiments performed at Studsvik, Sweden in the 1970s
  using light-water-moderated lattices with uranium rods (1.35% U-235 UO2) and
  mixed-oxide (MOX) rods. The KRITZ reactor used a pressure vessel large enough for
  multiple full-size fuel assemblies. Criticality achieved at room temperature (~20 C)
  and elevated temperature (~250 C) by controlling boron content and water level.
- **Specification URL:** https://www.oecd-nea.org/upload/docs/application/pdf/2019-12/nsc-doc2005-24.pdf
  (OECD/NEA benchmark report) and ICSBEP LEU-COMP-THERM-104.
- **Reference results:** Critical boron concentrations, critical water levels, axial
  buckling measurements, relative rod powers for selected fuel rods. Temperature-dependent
  k uncertainties reduced from 195 pcm to 40 pcm for large temperature changes.
- **Complexity:** ~500-1,500 surfaces for a single lattice configuration; multiple
  experimental cases available.
- **Notes:** Excellent for validating temperature-dependent reactivity (Doppler effect)
  modeling. The 2019 IRPhEP Handbook edition includes 37 accepted benchmark measurements
  for KRITZ-1-Mk. Freely available OECD/NEA specification PDF.

---

## 8. VVER-1000 Mock-up Critical Experiments

- **Type:** Critical experiment (hexagonal lattice PWR mock-up)
- **Description:** Mock-up of a VVER-1000 reactor with 32 dismountable fuel assemblies in
  a hexagonal lattice with 23.6 cm pitch, plus full-scale component simulators of
  baffle, barrel, displacer, pressure vessel, and biological shielding. Provides
  hexagonal fuel geometry distinct from Western square-lattice PWR designs.
- **Specification URL:** https://scholarworks.unist.ac.kr/bitstream/201301/48692/2/1_s2.0_S1738573320301819_main%20(1).pdf.pdf
  (open access validation paper)
- **Reference results:** Six experimental critical configurations, four pin-by-pin power
  maps. Code comparison between MCS and MCNP6 using ENDF/B-VII.1 and ENDF/B-VIII.0.
- **Complexity:** ~2,000-5,000 surfaces for the 32-assembly mock-up; a single hexagonal
  assembly is ~300-800 surfaces.
- **Notes:** Important for validating hexagonal geometry handling in OpenMC. The X2
  VVER-1000 benchmark (see doi:10.1016/j.anucene.2020.107698) provides additional
  fresh HZP core data with Serpent 2 reference Monte Carlo solutions. VVER geometry
  exercises OpenMC's hexagonal lattice capabilities.

---

## 9. SNAP-10A/2 Space Reactor Critical Experiments

- **Type:** Critical experiment (compact space reactor)
- **Description:** Criticality benchmarks for SNAP-10A reactor cores from the SCA-4B
  experimental program at Atomics International. The core is a cylindrical stainless
  steel vessel containing 37 closely spaced fuel elements (1.250 in diameter, 12.450 in
  long) on a triangular pitch of 1.260 in. Fuel is highly enriched uranium-zirconium
  hydride (UZrH). Various water immersion and reflection conditions tested.
- **Specification URL:** https://info.ornl.gov/sites/publications/Files/Pub57468.pdf
  (ORNL report, freely available)
- **Reference results:** Experimental k-eff for multiple configurations; OpenMC validation
  using ENDF/B-VII.0 and ENDF/B-VIII.0 cross section libraries.
- **Complexity:** ~100-300 surfaces; compact 37-element core is relatively simple.
- **Notes:** Already validated with OpenMC. Unique compact reactor type (space power).
  The uranium-zirconium hydride fuel and compact geometry provide a distinct validation
  case from water-moderated systems.

---

## 10. IAEA 10 MW MTR Research Reactor Benchmark

- **Type:** Research reactor benchmark (pool-type MTR)
- **Description:** An idealized pool-type 10 MW material testing reactor used as a
  standard benchmark for comparing computational methods. Specifications cover both
  HEU (93% U-235) and LEU (20% U-235) core configurations, supporting reactor
  conversion studies. The reactor uses plate-type fuel elements in a compact core
  arrangement with water moderator/coolant and beryllium/graphite reflectors.
- **Specification URL:** https://www-pub.iaea.org/MTCD/Publications/PDF/te_643v1_prn.pdf
  (IAEA-TECDOC-643, freely available)
- **Reference results:** k-eff, power distributions, thermal/fast flux distributions for
  six specified benchmark cases (BOL, EOL, HEU, LEU combinations). Burnup-dependent
  Monte Carlo calculations available.
- **Complexity:** ~500-1,500 surfaces; plate-type fuel elements with moderator channels.
- **Notes:** The IAEA MTR benchmark is one of the most widely reproduced reactor
  benchmarks in the world. Plate-type fuel geometry is distinct from pin-type
  assemblies. Multiple published validations exist across many codes.

---

## 11. OECD/NEA Burnup Credit Benchmark Phase I-B

- **Type:** PWR pin cell depletion benchmark
- **Description:** An infinite lattice of PWR fuel rods based on a Combustion Engineering
  14x14 assembly design. The benchmark specifies a simple pin cell with UO2 fuel,
  Zircaloy cladding, and water moderator. Actual pin dimensions are used with a modified
  pitch to match the assembly fuel-to-moderator ratio. Consists of 13 eigenvalue
  calculation cases investigating burnup, cooling time, and nuclide combinations.
- **Specification URL:** https://www.oecd-nea.org/science/wpncs/Publications/BUC/NSCDOC(96)06-buc-1b.pdf
  (OECD/NEA report NEA/NSC/DOC(96)-06, freely available)
- **Reference results:** k-inf for 13 cases with multiple burnup levels and cooling times;
  detailed isotopic compositions. Multi-code comparison results from international
  participants.
- **Complexity:** ~10-20 surfaces per pin cell (very simple geometry); the interest is in
  depletion accuracy rather than geometric complexity.
- **Notes:** Excellent for validating OpenMC's depletion capability. The Phase II-D
  benchmark extends this to include control rod insertion effects and additional burnup
  levels (30 and 45 GWd/tU). Simple geometry makes it easy to implement but physically
  meaningful validation of burnup/depletion.

---

## 12. ASTRA Critical Facility (Pebble Bed HTGR)

- **Type:** Critical facility experiment (pebble bed)
- **Description:** The ASTRA critical facility at Kurchatov Institute, Russia, is a
  standard benchmark for neutronics modeling of pebble bed HTGRs. Uses graphite pebbles
  with TRISO coated fuel particles in a cylindrical geometry. The benchmark specifies
  geometry, material compositions, and layout of spheres. The goal is to determine the
  critical pebble bed height.
- **Specification URL:** https://www.sciencedirect.com/science/article/abs/pii/S0306454921001675
  (OpenMC implementation paper)
- **Reference results:** Critical height measurements; OpenMC k-eff with standard
  deviations < 10 pcm. ENDF/B-VII.0 and VII.1 library comparisons; ENDF/B-VII.1
  outperforms VII.0 for criticality prediction.
- **Complexity:** ~500-1,500 surfaces; pebble packing and TRISO particle modeling
  required. Uses 2 million particles/batch, 150 batches (50 inactive).
- **Notes:** Already modeled in OpenMC with multiple geometry representations (core-hom,
  ref-hom, CR-hom). Complements HTR-10 as a second pebble bed benchmark. Well-suited
  for testing double heterogeneity treatment.

---

## 13. OECD/NEA SFR-UAM Sodium Fast Reactor Benchmark

- **Type:** Sodium-cooled fast reactor (pin cell to full core)
- **Description:** OECD/NEA benchmark for uncertainty analysis in sodium-cooled fast
  reactors. Two SFR core designs: a 3600 MWth oxide core and a 1000 MWth metallic core.
  For each core size, three fuel types are specified: oxide, carbide, and metal.
  Exercises range from pin cell to sub-assembly to full core under steady-state BOL and
  EOEC conditions.
- **Specification URL:** https://www.oecd-nea.org/jcms/pl_20438 (benchmark page;
  specifications accessible through OECD/NEA membership)
- **Reference results:** Doppler and void coefficients, dynamic characteristics, k-eff.
  Multi-code comparison results from international participants.
- **Complexity:** Pin cell is ~10-20 surfaces; single hexagonal sub-assembly is
  ~200-500 surfaces; the full core is much larger.
- **Notes:** Fast spectrum benchmark provides validation distinct from thermal reactor
  cases. The pin cell and single assembly exercises are tractable. Related benchmarks
  include EBR-II (SHRT-17, IAEA-TECDOC-1819) and FFTF (OpenMC validation in
  doi:10.1016/j.anucene.2023.110050). Access to full specifications may require
  OECD/NEA membership.

---

## 14. OECD/NEA LWR-UAM Benchmark Phase I (Neutronics)

- **Type:** LWR pin cell and assembly uncertainty benchmark
- **Description:** OECD/NEA benchmark for uncertainty analysis in best-estimate modeling
  of light water reactors. Phase I focuses on neutronics at three scales: pin cell,
  fuel lattice, and reactor core. Three reactor types are covered: PWR, BWR, and VVER.
  Each type has fully specified pin cell and assembly geometries with material compositions
  and operating conditions.
- **Specification URL:** https://www.oecd-nea.org/jcms/pl_80569 (OECD/NEA publication,
  Volume I: Neutronics Phase)
- **Reference results:** Multi-code comparison of k-inf, reaction rates, power
  distributions, and uncertainty propagation from nuclear data to integral parameters.
- **Complexity:** Pin cell: ~10-20 surfaces; assembly: ~200-500 surfaces; core: much larger.
  Each scale provides a distinct benchmark.
- **Notes:** Covers three distinct reactor types (PWR, BWR, VVER) in one benchmark
  framework. The pin cell and assembly exercises are immediately tractable. Provides
  both reference solutions and uncertainty quantification data. Access to full
  specifications may require OECD/NEA membership, but many participant papers publish
  the specifications openly.

---

## 15. OECD/NEA Hoogenboom-Martin Monte Carlo Performance Benchmark

- **Type:** Full-size PWR core performance benchmark
- **Description:** A simplified but full-size PWR model with 241 fuel assemblies, each a
  17x17 array of pins, with each pin subdivided into 100 axial nodes (over 7 million
  tally regions). Designed for benchmarking Monte Carlo code performance on detailed
  power density calculations. The geometry is a realistic PWR core but with simplified
  material definitions.
- **Specification URL:** https://www.researchgate.net/publication/242179791 (original paper);
  specifications also available in the MIT-CRPG benchmarks repository at
  https://github.com/mit-crpg/benchmarks (mc-performance directory includes OpenMC model).
- **Reference results:** Power density distributions with statistical uncertainties;
  40 billion neutron-history calculations achieving < 2% standard deviation in 43% of
  tally regions.
- **Complexity:** Full core is very large (~50,000+ surfaces), but a single assembly or
  mini-core subset is tractable (~300-1,000 surfaces). The OpenMC model already exists
  in the MIT-CRPG repository.
- **Notes:** Organized under OECD/NEA Data Bank auspices. OpenMC model already available.
  Best used as a scaling/performance benchmark rather than for physics validation, since
  the geometry is simplified. A single-assembly extraction provides a tractable physics
  problem.

---

## Summary Table

| # | Benchmark | Type | Geometry | Approx. Surfaces | OpenMC Model Exists? |
|---|-----------|------|----------|-------------------|---------------------|
| 1 | VERA Problems 1-3 | PWR pin/assembly | Square 17x17 | 10-500 | Yes |
| 2 | BEAVRS | PWR full core | Square 17x17 | 500-50,000+ | Yes |
| 3 | NuScale-Like SMR | SMR core | Square 17x17 | 300-10,000 | Yes |
| 4 | VENUS-2 MOX | Critical expt. | Square lattice, UO2+MOX | 1,000-3,000 | Yes |
| 5 | MSRE | Molten salt reactor | Cylindrical, graphite stringers | 500-2,000 | Yes |
| 6 | HTR-10 | Pebble bed HTGR | Cylindrical, TRISO pebbles | 500-2,000 | Partial |
| 7 | KRITZ-2 | Critical expt. | Square lattice, UO2+MOX | 500-1,500 | No |
| 8 | VVER-1000 Mock-up | Critical expt. | Hexagonal lattice | 300-5,000 | No |
| 9 | SNAP-10A/2 | Space reactor expt. | Cylindrical, triangular pitch | 100-300 | Yes |
| 10 | IAEA 10 MW MTR | Research reactor | Plate-type fuel | 500-1,500 | No |
| 11 | OECD/NEA BUC Phase I-B | PWR pin cell depletion | Single pin cell | 10-20 | No |
| 12 | ASTRA Pebble Bed | Critical facility | Cylindrical, pebble bed | 500-1,500 | Yes |
| 13 | SFR-UAM | Sodium fast reactor | Hexagonal, pin/assembly | 10-500 | No |
| 14 | LWR-UAM Phase I | LWR pin/assembly/core | Square + hexagonal | 10-500 | No |
| 15 | Hoogenboom-Martin | PWR full core | Square 17x17 | 300-50,000+ | Yes |

---

## Recommended Implementation Priority

**Tier 1 -- Immediate (existing OpenMC models, freely available specs):**
1. VERA Problem 2 (single 17x17 assembly)
2. NuScale-Like SMR (single assembly or small core)
3. SNAP-10A/2 (compact, simple geometry)
4. VENUS-2 MOX (mixed fuel critical experiment)

**Tier 2 -- High value, moderate effort:**
5. MSRE (unique reactor type, OpenMC model exists)
6. BEAVRS single assembly extraction
7. IAEA 10 MW MTR (plate fuel, widely published specs)
8. HTR-10 (TRISO/pebble bed, advanced reactor)

**Tier 3 -- Valuable but may need specification access:**
9. KRITZ-2 (temperature-dependent experiments)
10. VVER-1000 mock-up (hexagonal geometry validation)
11. OECD/NEA BUC Phase I-B (depletion validation)
12. ASTRA critical facility (pebble bed complement to HTR-10)
13. SFR-UAM pin cell/assembly (fast spectrum)
14. LWR-UAM pin cell/assembly (uncertainty quantification)
15. Hoogenboom-Martin single assembly (performance benchmark)

---

## Key Source Repositories

- **MIT-CRPG Benchmarks:** https://github.com/mit-crpg/benchmarks
- **BEAVRS:** https://github.com/mit-crpg/BEAVRS
- **ECP ExaSMR:** https://github.com/mit-crpg/ecp-benchmarks
- **ICSBEP Handbook:** https://www.oecd-nea.org/jcms/pl_20291
- **IRPhE Handbook:** https://www.oecd-nea.org/jcms/pl_20279
- **OECD/NEA Reactor Physics Benchmarks:** https://www.oecd-nea.org/science/projects/benchmarks.html
- **INL Virtual Test Bed:** https://mooseframework.inl.gov/virtual_test_bed/
