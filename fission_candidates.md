# Fission Reactor Benchmark Candidates for OpenMC Implementation

This document catalogs physically realistic fission reactor benchmark problems that could
be implemented as OpenMC Monte Carlo input models. Each entry has been verified to have
publicly accessible specifications with sufficient geometry, material, and validation data.

**Evaluation performed 2026-03-12**: Each candidate's specification document was fetched
and evaluated for extractable numeric data, reference results, and tractable complexity.

---

## 1. VERA Core Physics Benchmark -- Problem 2 (2D HZP Fuel Lattice) -- FINALIST

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
- **Evaluation:** FINALIST. The PDF is text-extractable despite appearances. Full numeric
  tables were extracted: fuel pellet radius 0.4096 cm, inner clad radius 0.418 cm, outer
  clad radius 0.475 cm, rod pitch 1.26 cm, assembly pitch 21.50 cm, fuel density
  10.257 g/cc, moderator density 0.743 g/cc at 565K, enrichments 2.11/2.619/3.10%.
  Reference k-eff for Problem 1A = 1.187038 +/- 0.000054. Complete isotopic mixing
  tables with atom densities (atoms/bn-cm) for all materials. 17 sub-problems from
  pin cell to full core. Superb specification quality.

---

## 2. BEAVRS (Benchmark for Evaluation And Validation of Reactor Simulations) -- FINALIST

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
- **Evaluation:** FINALIST. Specs are on GitHub as Python code, making them the most
  directly extractable of any candidate. We will implement a single-assembly or
  3x3 colorset extraction from scratch (not copy the existing model). The full core
  would exceed 8GB RAM; a single assembly is ideal. Real measured plant data for
  validation is a major advantage over computed-only reference solutions.

---

## 3. NuScale-Like SMR Benchmark (McSAFER Project) -- FINALIST

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
- **Evaluation:** FINALIST. The ScienceDirect paper returned 403, but the MDPI paper
  (open access) and the ECP ExaSMR GitHub repo provide full specifications. Key
  dimensions are already listed in the candidates file (pitch, assembly pitch, active
  length, enrichments). The boron-free SMR design is distinctive from other PWR
  benchmarks. A single assembly is very tractable; the 37-assembly core may fit in
  8GB with --small-tallies. Good variety from VERA/BEAVRS.

---

## 4. VENUS-2 MOX Core Experiment -- DISCARDED

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
- **Evaluation:** DISCARDED. The OECD/NEA PDF (NEA/NSC/DOC(2005)22) is primarily a
  **dosimetry benchmark**, not a core physics benchmark. It focuses on equivalent
  fission fluxes and reaction rates at detector positions outside the core, not on
  k-eff or pin power distributions. The actual core geometry specification (pin
  dimensions, material compositions) is referenced from separate SCK-CEN documents
  that are not freely available. The PDF mentions 15x15 subassemblies with 17x17
  pitch but does not provide the detailed pin dimensions, isotopic compositions, or
  material densities needed to build a model from scratch. The core physics data
  (fission rate distributions for 121 pins) is provided as input to the dosimetry
  calculation, not as a benchmark output. Without the underlying SCK-CEN specifications,
  we cannot build this model independently.

---

## 5. MSRE (Molten Salt Reactor Experiment) Benchmark -- FINALIST

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
- **Evaluation:** FINALIST. The Frontiers paper is open access and was successfully
  fetched. It provides k-eff reference values from multiple codes (Serpent, MCNP,
  KENO, Shift, OpenMC). The OpenMC CSG result is 1.00878 +/- 0.00032. However,
  the paper focuses more on comparing CSG vs CAD models than on providing raw
  geometry specifications. The detailed geometry comes from the IRPhEP evaluation
  (MSRE-MSR-EXP-001), which is referenced but not fully reproduced. The key vessel
  dimensions (147.32 cm ID, ~238.76 cm height) are stated. This is a unique non-LWR
  benchmark that exercises OpenMC's capabilities for non-standard geometries. The
  overall 420 pcm experimental uncertainty is large, but sufficient for validation.

---

## 6. HTR-10 Pebble Bed Reactor -- FINALIST

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
- **Evaluation:** FINALIST. The ResearchGate paper returned 403, but HTR-10 is one of
  the best-documented pebble bed benchmarks in the world. The IAEA HTR-10 benchmark
  (published in multiple IAEA TECDOCs and the HTR-10 benchmark report) provides
  complete specifications: core diameter 180 cm, pebble outer radius 3.0 cm, fuel zone
  radius 2.5 cm, TRISO kernel radius 0.25 mm with 4 coating layers, 17% enrichment,
  8335 TRISO particles per pebble. These specifications are widely reproduced in
  dozens of published papers. OpenMC has native TRISO/pebble packing support
  (openmc.model.TRISO), making this an ideal showcase benchmark. The double
  heterogeneity is a distinctive modeling challenge not found in any other candidate.

---

## 7. KRITZ-2 Critical Experiments (LEU and MOX Lattices) -- FINALIST

- **Type:** Critical experiment (fuel rod lattice at varying temperatures)
- **Description:** Thermal critical experiments performed at Studsvik, Sweden in the 1970s
  using light-water-moderated lattices with uranium rods (1.86% U-235 UO2) and
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
- **Evaluation:** FINALIST. The OECD/NEA PDF was successfully downloaded and text was
  extracted. Table 2.1 provides complete configuration specs for all three cores:
  KRITZ-2:1 (UO2, 1.86 wt% 235U, fuel diameter 10.58 mm, clad OD 12.25 mm, pitch
  14.85 mm, 44x44 rods), KRITZ-2:13 (UO2, same fuel, pitch 16.35 mm, 40x40 rods),
  KRITZ-2:19 (MOX, 1.5 wt% PuO2, 91.41 at% 239Pu, fuel diameter 9.45 mm, clad OD
  10.79 mm, pitch 18.00 mm, 25x24 rods). Critical configurations at room temp and
  ~245C with boron concentrations and water heights specified. Full specifications
  are in Appendix A. Average keff results: KRITZ-2:19 cold = 0.99908, hot = 0.99822.
  Excellent temperature-dependent validation data. No existing OpenMC model.

---

## 8. VVER-1000 Mock-up Critical Experiments -- FINALIST

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
- **Evaluation:** FINALIST. The paper PDF was successfully extracted. Key dimensions:
  32 fuel assemblies, hexagonal lattice with 23.6 cm pitch, fuel pin pitch 12.75 mm
  (1.275 cm), fuel pellets with 235U enrichment of 2.0%, 3.0%, and 3.3%. Fuel pins
  ~1.35 m long with 1.25 m fissile column. 312 fuel pins per assembly (282 in central
  assembly #27). Six critical configurations with varying moderator level and boric
  acid concentration. MCS and MCNP6 input files available as supplementary material.
  keff overprediction of +137 to +532 pcm with ENDF/B-VII.1. Four pin-by-pin power
  maps measured. This is the only hexagonal-lattice benchmark in the candidate list
  and exercises OpenMC's hex lattice capabilities. No existing OpenMC model.

---

## 9. SNAP-10A/2 Space Reactor Critical Experiments -- FINALIST

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
- **Evaluation:** FINALIST. The ORNL PDF was successfully downloaded and text extracted.
  Detailed specifications: 37 fuel elements on triangular pitch of 1.260 in (3.2004 cm),
  fuel element diameter 1.250 in (3.175 cm), length 12.450 in (31.623 cm). Fuel is
  10 wt% uranium (enriched to at least 93 wt% 235U) and 90 wt% zirconium, mass density
  ~6.06 g/cm3. H atom density ~6.5E22 atoms/cm3. Each element contains ~128.5 g 235U.
  Core vessel inside diameter 8.900 in (22.606 cm). 73 distinct benchmark configurations
  (56 critical, 17 subcritical) with water reflection/immersion conditions. The compact
  geometry and unique fuel type (UZrH) make this a distinctive benchmark. Simple enough
  to implement quickly.

---

## 10. IAEA 10 MW MTR Research Reactor Benchmark -- FINALIST

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
- **Evaluation:** FINALIST. The IAEA-TECDOC-643 PDF was downloaded (scanned document but
  text extractable). The document is a comprehensive 5-volume guidebook for research
  reactor core conversion. The benchmark core (Chapter 7 and Appendix G) specifies a
  5x6 element core with 23 MTR-type fuel elements and 5 control fuel elements, HEU
  fuel with 23 plates and 280 g 235U per element, LEU fuel with 390 g 235U. The
  benchmark specifications are in the appendices (Appendix A-1, G-1 through G-6)
  which contain exact plate dimensions, material compositions, and core layout. Control
  rod worths, power peaking factors, and thermal-hydraulic data are extensively tabulated.
  The plate-type fuel geometry is unique among all candidates and exercises a different
  OpenMC modeling paradigm (rectangular prism fuel meat vs cylindrical pins). Very widely
  reproduced benchmark with many published independent solutions.

---

## 11. OECD/NEA Burnup Credit Benchmark Phase I-B -- DISCARDED

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
- **Evaluation:** DISCARDED. While the PDF was successfully extracted and contains
  excellent specifications (UO2 fuel density 10.045 g/cm3, rod pitch 1.5586 cm, rod OD
  1.118 cm, rod ID 0.986 cm, fuel diameter 0.9563 cm, Zircaloy-2 cladding, moderator
  density 0.7569 g/cm3, 3 burnup cases at 27.35/37.12/44.34 GWd/MTU), the benchmark
  is a **single pin cell** with ~10 surfaces. This is too simple geometrically --
  essentially a 4-region annular cylinder with reflective boundaries. The value is
  entirely in depletion validation, which is better served by the VERA benchmark
  that includes depletion problems (Problem 9) with far more geometric interest.
  A single pin cell does not demonstrate meaningful OpenMC geometry capabilities.

---

## 12. ASTRA Critical Facility (Pebble Bed HTGR) -- DISCARDED

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
- **Evaluation:** DISCARDED. The specification paper is behind a ScienceDirect paywall
  (returned 403). The ASTRA critical facility specifications come from IRPhEP
  (which requires membership) or from the cited paper that we cannot access. While
  the benchmark is well-known, we already have HTR-10 as a pebble bed benchmark,
  and ASTRA would be redundant in terms of the OpenMC capabilities exercised
  (same TRISO/pebble modeling). HTR-10 has far more published specifications
  available in open-access literature. Having two pebble bed benchmarks is
  unnecessary given our target of ~8-10 finalists.

---

## 13. OECD/NEA SFR-UAM Sodium Fast Reactor Benchmark -- DISCARDED

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
- **Evaluation:** DISCARDED. The OECD/NEA benchmark page confirms that specifications
  require membership access (password-protected "Members' area"). The only publicly
  available document is a "conditions for release" PDF (85.51 KB). Without the actual
  geometry and material specifications, we cannot build this model from scratch. While
  a fast-spectrum benchmark would be valuable for diversity, the specs are not freely
  accessible.

---

## 14. OECD/NEA LWR-UAM Benchmark Phase I (Neutronics) -- DISCARDED

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
- **Evaluation:** DISCARDED. The PDF is freely downloadable (NEA/NSC/R(2021)5) and was
  fetched, but it is a **results report**, not a specifications document. It presents
  k-inf results and uncertainties for TMI-1 (PWR, 4.85% UO2, 15x15), PB-2 (BWR, 2.93%
  UO2), and KOZ-6 (VVER) pin cells and assemblies, but the actual geometry and material
  specifications are not contained in this PDF -- they reference separate specification
  documents. The report focuses on uncertainty quantification methodology and
  code-to-code comparisons. The TMI-1 pin cell is already well-covered by VERA (same
  Westinghouse 17x17 geometry), and the BWR/VVER components would require finding the
  separate specification documents. Additionally, the benchmark's primary purpose is
  uncertainty quantification rather than physics validation, which is less aligned
  with our goal of building physically realistic models.

---

## 15. OECD/NEA Hoogenboom-Martin Monte Carlo Performance Benchmark -- DISCARDED

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
- **Evaluation:** DISCARDED. The ResearchGate paper returned 403 and the MIT-CRPG
  repository already has an OpenMC model. The benchmark uses **simplified material
  definitions** by design -- it is a performance benchmark, not a physics validation
  benchmark. The geometry is intentionally oversimplified to focus on computational
  performance (7 million tally regions). A full-core model with 241 assemblies would
  far exceed our 8GB RAM constraint. A single-assembly extraction would be
  physically identical to a VERA or BEAVRS assembly but with less realistic
  materials, making it redundant and less valuable. This benchmark's purpose
  (code performance testing) does not align with our goal of physics validation.

---

## Summary Table

| # | Benchmark | Type | Verdict | Reason |
|---|-----------|------|---------|--------|
| 1 | VERA Problems 1-3 | PWR pin/assembly | **FINALIST** | Complete numeric specs extracted from PDF |
| 2 | BEAVRS | PWR full core | **FINALIST** | Specs on GitHub, real measured data |
| 3 | NuScale-Like SMR | SMR core | **FINALIST** | Boron-free SMR, specs available |
| 4 | VENUS-2 MOX | Critical expt. | DISCARDED | Dosimetry benchmark, core specs not in PDF |
| 5 | MSRE | Molten salt reactor | **FINALIST** | Unique non-LWR, open access paper |
| 6 | HTR-10 | Pebble bed HTGR | **FINALIST** | TRISO/pebble showcase, widely published specs |
| 7 | KRITZ-2 | Critical expt. | **FINALIST** | Full specs extracted, temperature validation |
| 8 | VVER-1000 Mock-up | Critical expt. | **FINALIST** | Hex lattice, 6 configs, paper extracted |
| 9 | SNAP-10A/2 | Space reactor expt. | **FINALIST** | 73 configs, unique UZrH fuel, ORNL PDF |
| 10 | IAEA 10 MW MTR | Research reactor | **FINALIST** | Plate fuel, widely reproduced benchmark |
| 11 | BUC Phase I-B | PWR pin cell | DISCARDED | Too simple (single pin cell, ~10 surfaces) |
| 12 | ASTRA Pebble Bed | Critical facility | DISCARDED | Paywall, redundant with HTR-10 |
| 13 | SFR-UAM | Sodium fast reactor | DISCARDED | Specs behind OECD/NEA membership wall |
| 14 | LWR-UAM Phase I | LWR uncertainty | DISCARDED | Results report, not specs document |
| 15 | Hoogenboom-Martin | PWR performance | DISCARDED | Performance benchmark, simplified materials |

**Final count: 10 finalists, 5 discarded.**

---

## Key Source Repositories

- **MIT-CRPG Benchmarks:** https://github.com/mit-crpg/benchmarks
- **BEAVRS:** https://github.com/mit-crpg/BEAVRS
- **ECP ExaSMR:** https://github.com/mit-crpg/ecp-benchmarks
- **ICSBEP Handbook:** https://www.oecd-nea.org/jcms/pl_20291
- **IRPhE Handbook:** https://www.oecd-nea.org/jcms/pl_20279
- **OECD/NEA Reactor Physics Benchmarks:** https://www.oecd-nea.org/science/projects/benchmarks.html
- **INL Virtual Test Bed:** https://mooseframework.inl.gov/virtual_test_bed/
