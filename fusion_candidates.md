# Fusion Device Benchmark Candidates for OpenMC Input Models

This document catalogs fusion neutronics benchmark problems that could be implemented
as OpenMC Monte Carlo input models. All candidates have detailed geometry/material
specifications and reference experimental or computational results for validation.

---

## 1. FNG-ITER Bulk Shield Mock-up

- **Type:** Blanket/shield mock-up experiment
- **Description:** A 94 cm thick mock-up of the ITER inboard shielding (first wall,
  blanket, vacuum vessel, and toroidal field coils) irradiated by 14 MeV neutrons at
  the ENEA Frascati Neutron Generator. The assembly consists of copper, stainless
  steel SS316, and Perspex (water-equivalent) sandwich layers, with a rear section of
  alternating Cu/SS plates simulating the magnet.
- **Specification URL:** https://www.oecd-nea.org/science/wprs/shielding/sinbad/FNG_BLKT/FNGBKT_A.HTM
- **Reference results available:**
  - Activation foil reaction rates at 14 depths (3.4--92.2 cm): Nb-93(n,2n),
    Al-27(n,alpha), Ni-58(n,p), Au-197(n,gamma)
  - Nuclear heating via TLD-300 dosimeters at 17 positions
  - MCNP calculated C/E ratios with FENDL-1 and EFF-3
- **Data files:** MCNP-4A geometry input, FORTRAN source routine, experimental
  reaction rate and nuclear heating tables
- **Complexity:** Moderate -- slab geometry with multiple material layers; ~100 cells
- **Notes:** One of the most widely used SINBAD fusion benchmarks. Already used for
  TRIPOLI-4 vs MCNP numerical benchmarking. Well-documented in SINBAD.

---

## 2. FNG-ITER Streaming Experiment

- **Type:** Neutron streaming/shielding experiment
- **Description:** A neutron streaming benchmark through an ITER-like shielding
  assembly with a 28 mm diameter duct channel, a detector cavity, and rear
  Cu/SS316 plates simulating tokamak coils. Two configurations tested: on-axis
  (source aligned with channel) and off-axis (5.3 cm lateral shift). Front cross
  section 100x100 cm, total thickness 94.26 cm.
- **Specification URL:** https://www.oecd-nea.org/science/wprs/shielding/sinbad/fng_str/fngstr-a.htm
- **Reference results available:**
  - Activation foil reaction rates (Nb-93, Al-27, Ni-58, Au-197) at positions
    in channel, cavity, and behind assembly
  - Nuclear heating via TLD-300 dosimeters
  - MCNP results with FENDL-1/2 and EFF-3 libraries
- **Data files:** Full 3D MCNP input files (mcnpfoil.inp at 73 KB, mcnp_nh.inp at
  87 KB), FORTRAN source routine, experimental tables
- **Complexity:** Moderate-high -- 3D geometry with streaming channel and cavity;
  tests weight-window variance reduction
- **Notes:** Already validated with OpenMC (published 2024 in Fusion Science and
  Technology). Excellent candidate since OpenMC results already exist for comparison.
  Ideal test for weight windows / random ray solver.

---

## 3. FNG HCPB Tritium Breeder Module Mock-up

- **Type:** Breeding blanket mock-up experiment
- **Description:** A mock-up of the Helium-Cooled Pebble Bed (HCPB) breeding blanket
  concept, consisting of a metallic beryllium assembly with two double-layer breeder
  sections of Li2CO3 powder in an AISI-303 stainless steel box (31x29x30.9 cm).
  Rear cassette of AISI-316 steel holds additional Li2CO3. Irradiated by 14 MeV
  neutrons at FNG.
- **Specification URL:** https://www.oecd-nea.org/science/wprs/shielding/sinbad/fng_hcpb/fnghcpb-a.htm
- **Reference results available:**
  - Tritium production rates from Li-6(n,t) and Li-7(n,t) at 16 positions
  - Neutron reaction rates (Au, Ni, Al, Nb foils) at 4 positions in beryllium
  - Nuclear heating via TLD-300 in breeder layers
  - MCNP-4C, DORT, TORT analysis results
- **Data files:** MCNP-4C model (mcnp-hcpb.i), D-T source routines for MCNP5/MCNPX,
  DORT/TORT input files, geometry figures, 32 total files
- **Complexity:** Moderate -- box geometry with layered internals; detailed material
  compositions
- **Notes:** Ranked as benchmark-quality experiment. Already converted to Serpent and
  OpenMC inputs in published benchmarking studies. Key for tritium breeding validation.

---

## 4. FNG Copper Benchmark

- **Type:** Material shielding experiment
- **Description:** A pure copper assembly (60x60x69.9 cm, seven 10-cm plates)
  irradiated by 14 MeV neutrons at FNG (2014-2015). Designed to validate copper
  nuclear cross-section data relevant for ITER design, as copper is a key structural
  material in ITER.
- **Specification URL:** https://www.sciencedirect.com/science/article/abs/pii/S0920379618308184
- **Reference results available:**
  - Reaction rates for 7 foil types (Nb-93, Al-27, Ni-58, In-115, Au-197(n,2n),
    W-186(n,gamma), Au-197(n,gamma)) vs penetration depth
  - Neutron flux spectra (NE213 scintillator)
  - Dose rates (TLD)
  - MCNP analysis with multiple nuclear data libraries
- **Data files:** Full and simplified MCNP models available; entered into SINBAD
- **Complexity:** Low-moderate -- simple slab geometry of pure material
- **Notes:** Already converted to OpenMC and Serpent inputs in published benchmarking
  studies. One of the newer SINBAD entries with modern uncertainty quantification.

---

## 5. FNG Tungsten Experiment

- **Type:** Material shielding experiment
- **Description:** Measurements of neutron reaction rates and gamma heating in a
  tungsten assembly (DENSIMET-176/180 alloy) irradiated by 14 MeV neutrons at FNG.
  Tungsten is the leading plasma-facing material candidate for fusion reactors.
- **Specification URL:** https://www.oecd-nea.org/science/wprs/shielding/sinbad/fng_w/fngw-a.htm
- **Reference results available:**
  - Neutron reaction rates vs depth
  - Absorbed dose via TLD
  - MCNP C/E ratios
- **Data files:** Available in SINBAD (generated 2004, updated 2006)
- **Complexity:** Low-moderate -- slab/block geometry
- **Notes:** Tungsten cross-section validation is critical for ITER/DEMO divertor
  design. Also available as FNG/TUD variant with spectral measurements.

---

## 6. FNG-ITER Dose Rate Experiment

- **Type:** Shutdown dose rate experiment
- **Description:** A SS316/Perspex assembly (100x100 cm front, 71.83 cm thick) with a
  cavity and streaming channel, designed to validate shutdown dose rate calculations
  for ITER vacuum vessel conditions. Irradiated for 18 hours by 14 MeV FNG neutrons
  (May 2000), with dose rate measurements from 0.5 hours to 3+ months cooling time.
- **Specification URL:** https://www.oecd-nea.org/science/wprs/shielding/sinbad/fng_dose/fngdos-a.htm
- **Reference results available:**
  - Continuous dose rate decay curves (Geiger-Mueller detector)
  - Integrated dose at 4 cooling times (TLD-300, GR-200A)
  - Activation foil rates (Ni-58(n,p), Ni-58(n,2n))
  - R2S analysis with MCNP-4C + FISPACT
- **Data files:** 3D MCNP neutron transport model (mcnp_n.inp), decay gamma transport
  model (mcnp_g.inp), FISPACT inputs, FORTRAN source routine
- **Complexity:** Moderate -- 3D geometry with cavity; two-step R2S analysis needed
  for shutdown dose
- **Notes:** Tests activation/shutdown dose workflow (requires coupling OpenMC with an
  inventory code like FISPACT or OpenMC depletion). Already benchmarked with OpenMC.

---

## 7. FNS Clean Benchmark: Tungsten, Vanadium, and Beryllium Assemblies

- **Type:** Material neutronics experiments
- **Description:** Three clean benchmark experiments performed at the JAEA Fusion
  Neutronics Source (FNS) facility: a tungsten cylindrical assembly (diameter 629 mm,
  height 507 mm), a vanadium cube (25.4 cm side), and a beryllium assembly, all
  irradiated with 14 MeV D-T neutrons. Measurements at multiple depths.
- **Specification URLs:**
  - Tungsten: https://www.oecd-nea.org/science/wprs/shielding/sinbad/fns_w/fnsw-abs.htm
  - Vanadium: https://www.oecd-nea.org/science/wprs/shielding/sinbad/fns_v/fnsv-abs.htm
  - General FNS: https://www.oecd-nea.org/science/wprs/shielding/sinbad/fns/fns.htm
- **Reference results available:**
  - Neutron spectra (NE213 + proton recoil counters)
  - Photon spectra
  - Activation foil reaction rates at 4 positions
  - Gamma heating rates (TLD)
  - Tritium production rates (beryllium case)
  - MCNP-4A reference calculations
- **Data files:** MCNP input files, tabulated spectra, reference JAERI-Data reports
- **Complexity:** Low -- simple cylindrical/cubic geometry, single material
- **Notes:** Already fully validated with OpenMC v0.14.1-dev (published 2024). These
  are excellent first benchmarks due to simple geometry and comprehensive data.
  Available in the openmc_fusion_benchmarks GitHub repository.

---

## 8. FNS Dogleg Duct Streaming Experiment

- **Type:** Neutron streaming experiment
- **Description:** An iron slab assembly (1700x1400x1800 mm) with a doubly bent
  (dogleg) duct of 300x300 mm cross-section. Duct legs measure 1150, 600, and
  650 mm. D-T neutron source at FNS facility (~4E12 n/s). Simulates neutron
  streaming through diagnostic/heating channels in fusion reactor shielding.
- **Specification URL:** https://www.oecd-nea.org/science/wprs/shielding/sinbad/fns_duct/fnsstr-a.htm
- **Reference results available:**
  - Neutron spectra above 2 MeV (NE213 scintillator) at 4 positions
  - Activation dosimetry: Nb-93(n,2n), In-115(n,n'), Au-197(n,gamma)
  - MCNP-4B/4C calculations with FENDL/2 and JENDL-3.3
- **Data files:** Two MCNP input files (FENDL/2 and JENDL-3.3 libraries), measured
  and calculated data in Excel format
- **Complexity:** Moderate -- iron slab with bent duct geometry; streaming problem
  benefits from weight windows
- **Notes:** Important for validating variance reduction techniques. Being implemented
  in the openmc_fusion_benchmarks repository.

---

## 9. KANT Beryllium Spherical Shell Experiment

- **Type:** Material neutron multiplication experiment
- **Description:** The Karlsruhe Neutron Transmission (KANT) experiment measured
  neutron leakage spectra through concentric spherical beryllium shells of 5, 10,
  and 17 cm thickness (inner diameter 10 cm, outer up to 44 cm) with a central
  T(d,n) source. Measurements at 60 degrees emission angle.
- **Specification URL:** https://www.oecd-nea.org/science//wprs/shielding/sinbad/kant/fzk-be_a.htm
- **Reference results available:**
  - 178-group neutron leakage spectra (50 keV to 15 MeV) for 3 thicknesses
  - Five-group partial leakage multiplications
  - Multiple detector types: NE-213, proton recoil counters, TOF, Bonner spheres
  - 8% combined uncertainty (1-sigma)
- **Data files:** MCNP model input (kantmcnp.i), 178-group tabulated spectra
  (fzkbe.tb1-3), measured spectra spreadsheet (kant.xls), reference paper
- **Complexity:** Low -- spherical shell geometry, single material
- **Notes:** Critical for validating beryllium neutron multiplication data used in
  breeding blanket design. Analyzed with TRIPOLI-4, MCNP, and other codes.

---

## 10. OKTAVIAN Sphere Experiments (Iron, Nickel, Aluminium, Silicon, Tungsten, Manganese)

- **Type:** Material leakage spectra experiments
- **Description:** A series of sphere experiments performed at Osaka University's
  OKTAVIAN facility using a 14 MeV D-T neutron source at the center of various
  material spheres. The iron sphere has radius 50.32 cm, nickel sphere 16 cm
  diameter. Neutron leakage spectra measured by time-of-flight technique. Six
  different single-material sphere experiments are available.
- **Specification URLs:**
  - Iron: https://www.oecd-nea.org/science/wprs/shielding/sinbad/oktav_fe/okfe-abs.htm
  - Nickel: https://www.oecd-nea.org/science/wprs/shielding/sinbad/oktav_ni/okni-abs.htm
- **Reference results available:**
  - Neutron leakage spectra (30 keV to 15 MeV) via time-of-flight
  - Multiple materials available for systematic validation
- **Data files:** Available in SINBAD; implemented in the openmc_fusion_benchmarks
  GitHub repository (e.g., oktavian_al model)
- **Complexity:** Very low -- single sphere with central point source
- **Notes:** Simplest possible fusion benchmarks. Ideal as introductory/validation
  models. Already implemented in the openmc_fusion_benchmarks repository with
  postprocessing notebooks.

---

## 11. IPPE Iron Spherical Shell Transmission Experiment

- **Type:** Material transmission experiment
- **Description:** Neutron transmission measurements through five spherical iron
  shells of different sizes (wall thicknesses 2.5 to 28 cm, radii 4.5 to 30 cm)
  with a central 14.1 MeV D-T source. Time-of-flight detector at 6.8 m, 8 degrees
  from beam axis. Performed at IPPE, Russia (1989-1995).
- **Specification URL:** https://www.oecd-nea.org/science/wprs/shielding/sinbad/ippe-fe/ippe_fe-a.htm
- **Reference results available:**
  - Neutron leakage spectra (50 keV to 15 MeV) for 5 shell thicknesses
  - MCNP-4C reference calculations
  - Detailed material compositions from measured weights
- **Data files:** 5 MCNP-4C input files (one per sphere), experimental description
  document, reference report (EFF-DOC-747)
- **Complexity:** Very low -- spherical shell geometry, single material, 5 variants
- **Notes:** Provides a systematic thickness study allowing validation of deep
  penetration transport. Complementary to OKTAVIAN single-sphere experiments.

---

## 12. FNG WCLL (Water-Cooled Lithium Lead) Mock-up Experiment

- **Type:** Breeding blanket mock-up experiment
- **Description:** A mock-up of the Water-Cooled Lithium Lead (WCLL) breeding blanket
  concept for EU DEMO, consisting of LiPb bricks, EUROFER plates, and Perspex
  (substituting water). Irradiated by 14 MeV neutrons at FNG. Tritium production
  rate and detector reaction rates measured using Li2CO3 pellets and activation
  foils at positions up to ~55 cm depth.
- **Specification URL:** https://www.sciencedirect.com/science/article/abs/pii/S0920379620301484
- **Reference results available:**
  - Tritium production rates at multiple positions
  - Activation foil reaction rates vs depth
  - Sensitivity/uncertainty analysis with SUSD3D
  - Pre-analysis and post-analysis MCNP/deterministic calculations
- **Data files:** Detailed MCNP model; sensitivity/uncertainty data
- **Complexity:** Moderate -- multi-material layered assembly with realistic blanket
  materials (LiPb, EUROFER, water-equivalent)
- **Complexity:** Moderate -- multi-material layered geometry
- **Notes:** Most recent FNG blanket experiment (2019). Tests the WCLL concept which
  is one of two EU DEMO blanket candidates. Excellent for tritium breeding validation.

---

## 13. EU DEMO HCPB Breeding Blanket Sector Model

- **Type:** Computational tokamak sector benchmark
- **Description:** A 11.25-degree toroidal sector model of the EU DEMO reactor with
  Helium-Cooled Pebble Bed (HCPB) breeding blanket. Includes inboard and outboard
  blanket segments, vacuum vessel, and toroidal field coils. Uses reflecting boundary
  conditions for full-torus extrapolation. Based on EU DEMO BL2017 reference design.
- **Specification URL:** https://www.sciencedirect.com/science/article/abs/pii/S0920379620301319
- **Reference results available:**
  - Tritium breeding ratio (TBR)
  - Nuclear heating distributions
  - Neutron/photon flux spectra
  - Serpent-2 vs MCNP5 benchmark comparison
- **Data files:** Serpent-2 and MCNP models described in literature; geometry
  converted using csg2csg tool
- **Complexity:** High -- full tokamak sector with heterogeneous blanket modules,
  multiple material zones
- **Notes:** Not a physical experiment but a computational benchmark between codes.
  Represents a physically realistic reactor design. May require significant effort
  to build from scratch in OpenMC, but the csg2csg conversion tool could help.

---

## 14. TUD Iron Slab Experiment

- **Type:** Material shielding experiment
- **Description:** Spectral neutron and photon flux measurements from iron assemblies
  irradiated with 14 MeV neutrons at Technische Universitaet Dresden. Geometry
  chosen to be sensitive to shield penetration problems relevant to fusion reactor
  design. Measurements include neutron and gamma spectra at multiple depths.
- **Specification URL:** https://www.oecd-nea.org/science/wprs/shielding/sinbad/tud_fe/tufe-abs.htm
- **Reference results available:**
  - Energy spectra (neutron and photon) with statistical and systematic uncertainties
  - MCNP calculations with FENDL-1 and EFF-2 libraries
  - Tabulated energy, spectrum, and uncertainty data
- **Data files:** Available in SINBAD (generated 1998, updated 2010)
- **Complexity:** Low -- slab geometry, single material
- **Notes:** Widely used as one of three SINBAD fusion deep-penetration benchmarks for
  Monte Carlo code validation (alongside FNG bulk shield and IPPE iron).

---

## 15. Juelich Lithium Metal Blanket Experiment

- **Type:** Breeding blanket experiment
- **Description:** Blanket neutronics experiments performed at Institut fuer
  Reaktorentwicklung (IRE) in Juelich (1976-1984) using a 14 MeV D-T source.
  Multiple configurations tested: Li only, Be-Li (with beryllium neutron multiplier),
  Li-C (with graphite reflector), and Be-Li-C (both). Tritium production and
  neutron spectra measured with Li2CO3 samples, TLD detectors, and activation foils.
- **Specification URL:** https://www.oecd-nea.org/science/wprs/shielding/sinbad/juelich_li/Juelich_Li_met_a.htm
- **Reference results available:**
  - Space-dependent neutron spectra
  - Energy deposition profiles
  - Tritium production rates in multiple configurations
  - Effect of beryllium multiplier and graphite reflector quantified
- **Data files:** Available in SINBAD (generated/updated 2008)
- **Complexity:** Moderate -- cylindrical/annular blanket geometry with multiple
  material configurations
- **Notes:** Tests the fundamental physics of tritium breeding: lithium breeding,
  beryllium multiplication, and graphite reflection. Four configurations from a
  single facility provide systematic validation.

---

## 16. ORNL 14 MeV SS/Borated Polyethylene Slab Experiment

- **Type:** Shielding experiment
- **Description:** Integral experiments measuring neutron and gamma-ray energy spectra
  from the transport of ~14 MeV T(d,n) neutrons through laminated stainless steel
  and borated-polyethylene shield configurations at Oak Ridge National Laboratory
  (1979). NE-213 detector with pulse-shape discrimination used for spectral
  measurements.
- **Specification URL:** https://www.oecd-nea.org/science/wprs/shielding/sinbad/sinbadis.htm
  (entry in SINBAD index)
- **Reference results available:**
  - Neutron and gamma-ray energy spectra behind shield
  - Multiple shield configurations (laminated slabs)
- **Data files:** Available through SINBAD/RSICC
- **Complexity:** Low -- laminated slab geometry
- **Notes:** Historical US benchmark. Tests multi-material shielding with boron
  capture effects. Less commonly used than FNG/FNS experiments but provides
  independent US experimental data.

---

## Summary Table

| # | Name | Type | Complexity | OpenMC Prior Work | Experimental Data |
|---|------|------|------------|-------------------|-------------------|
| 1 | FNG-ITER Bulk Shield | Shield mock-up | Moderate | TRIPOLI/MCNP | Yes (SINBAD) |
| 2 | FNG-ITER Streaming | Streaming/shielding | Moderate-High | Yes (published) | Yes (SINBAD) |
| 3 | FNG HCPB TBM Mock-up | Breeding blanket | Moderate | Yes (published) | Yes (SINBAD) |
| 4 | FNG Copper | Material shielding | Low-Moderate | Yes (published) | Yes (SINBAD) |
| 5 | FNG Tungsten | Material shielding | Low-Moderate | No | Yes (SINBAD) |
| 6 | FNG-ITER Dose Rate | Shutdown dose | Moderate | Yes (published) | Yes (SINBAD) |
| 7 | FNS Clean (W/V/Be) | Material neutronics | Low | Yes (published) | Yes (SINBAD) |
| 8 | FNS Dogleg Duct | Streaming | Moderate | In progress | Yes (SINBAD) |
| 9 | KANT Be Shells | Neutron multiplication | Low | No | Yes (SINBAD) |
| 10 | OKTAVIAN Spheres | Leakage spectra | Very Low | Yes (repository) | Yes (SINBAD) |
| 11 | IPPE Iron Shells | Transmission | Very Low | No | Yes (SINBAD) |
| 12 | FNG WCLL Mock-up | Breeding blanket | Moderate | No | Yes |
| 13 | EU DEMO HCPB Sector | Tokamak sector | High | Partial | Computational |
| 14 | TUD Iron Slab | Material shielding | Low | No | Yes (SINBAD) |
| 15 | Juelich Li Blanket | Breeding blanket | Moderate | No | Yes (SINBAD) |
| 16 | ORNL SS/BPoly Slab | Shielding | Low | No | Yes (SINBAD) |

---

## Recommended Implementation Priority

**Phase 1 -- Simple geometries (weeks):**
- OKTAVIAN spheres (#10) -- simplest possible, already in openmc_fusion_benchmarks
- IPPE iron shells (#11) -- 5 spherical variants, MCNP inputs available
- KANT beryllium shells (#9) -- spherical, tests neutron multiplication
- FNS clean benchmarks (#7) -- already validated with OpenMC

**Phase 2 -- Moderate complexity (weeks to months):**
- FNG copper (#4) -- already benchmarked with OpenMC
- FNG tungsten (#5) -- similar to copper, important for ITER
- TUD iron slab (#14) -- slab geometry, deep penetration
- FNG HCPB TBM (#3) -- key breeding blanket benchmark

**Phase 3 -- Complex/streaming problems (months):**
- FNG-ITER bulk shield (#1) -- multi-material ITER mock-up
- FNG-ITER streaming (#2) -- already validated with OpenMC, tests weight windows
- FNS dogleg duct (#8) -- streaming through bent duct
- FNG-ITER dose rate (#6) -- requires R2S coupling for shutdown dose

**Phase 4 -- Blanket systems and reactor models (months):**
- FNG WCLL mock-up (#12) -- latest blanket experiment
- Juelich Li blanket (#15) -- multiple blanket configurations
- EU DEMO HCPB sector (#13) -- full reactor sector model

---

## Key Resources

- **SINBAD Database (NEA):** https://www.oecd-nea.org/science/wprs/shielding/sinbad/
- **SINBAD Index:** https://www.oecd-nea.org/science/wprs/shielding/sinbad/sinbadis.htm
- **openmc_fusion_benchmarks (GitHub):** https://github.com/eepeterson/openmc_fusion_benchmarks
- **OpenMC FNS validation paper:** https://www.tandfonline.com/doi/full/10.1080/15361055.2024.2323747
- **OpenMC FNG-Streaming validation:** https://www.tandfonline.com/doi/full/10.1080/15361055.2024.2400762
- **Emergent codes benchmarking paper:** https://www.sciencedirect.com/science/article/pii/S0920379622001934
- **SINBAD overview paper:** https://www.sciencedirect.com/science/article/pii/S0306454921001304
- **Serpent-2 vs MCNP DEMO HCPB:** https://www.sciencedirect.com/science/article/abs/pii/S0920379620301319
- **RSICC SINBAD distribution:** https://rsicc.ornl.gov/codes/dlc/dlc2/dlc-237.html
