# Fusion Device Benchmark Candidates for OpenMC Input Models

This document catalogs fusion neutronics benchmark problems that could be implemented
as OpenMC Monte Carlo input models. All candidates have detailed geometry/material
specifications and reference experimental or computational results for validation.

Each candidate has been evaluated against specification documents and tagged as
**FINALIST** or **DISCARDED** based on: availability of geometry specs, material
compositions, MCNP inputs, experimental reference data, and tractable complexity.

---

## 1. FNG-ITER Bulk Shield Mock-up -- FINALIST

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
- **Evaluation:** FINALIST. SINBAD page provides 14 foil positions, 17 TLD positions,
  source-to-mockup distance (5.3 cm), 60-degree source aperture, and MCNP-4A input.
  Material compositions (Cu, SS316, Perspex) are standard. Excellent experimental
  data with multiple reaction types. One of the most-cited fusion shielding benchmarks.

---

## 2. FNG-ITER Streaming Experiment -- FINALIST

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
- **Evaluation:** FINALIST. Exceptionally well-documented: 6 MCNP input files
  (mcnpfoil.inp, mcnp_nh.inp, mcnp_hss.inp, mcnp_hcu.inp, mcnp_tld.inp),
  FORTRAN source routine. Channel diameter 28 mm, cavity 52x148x48 mm, rear coil
  block 47x47 cm. Foil positions at 23 depths. Already validated with OpenMC.
  Best candidate for weight window / variance reduction testing.

---

## 3. FNG HCPB Tritium Breeder Module Mock-up -- FINALIST

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
- **Evaluation:** FINALIST. SINBAD page provides exact box dimensions (31x29x30.9 cm
  external, 0.5 cm wall), rear box (31x14.8x30.9 cm), breeder layer thickness
  (1.2 cm), Be density (1.85 g/cm3), Li2CO3 densities (1.123 and 0.9413 g/cm3),
  Li isotopics (7.5 at% Li-6, 92.5 at% Li-7). MCNP-4C model (84 KB) available.
  16 tritium production measurement positions. Key for breeding blanket validation.

---

## 4. FNG Copper Benchmark -- FINALIST

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
- **Evaluation:** FINALIST. Although the primary paper is behind the ScienceDirect
  paywall (403 error on fetch), this experiment is entered into SINBAD with full
  and simplified MCNP models. Already converted to OpenMC in published studies,
  so exact geometry/material data is available through the openmc_fusion_benchmarks
  repository. Simple 7-plate pure copper slab -- trivial geometry with 7 foil types.

---

## 5. FNG Tungsten Experiment -- FINALIST

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
- **Evaluation:** FINALIST. SINBAD page provides excellent detail: DENSIMET-176
  (92.3% W, 2.6% Fe, 4.2% Ni) and DENSIMET-180 (95.0% W, 1.6% Fe, 3.4% Ni)
  compositions, block dimensions (42-47x46.85x49 cm), source distance 5.3 cm,
  TLD-300 dosimeter specs. Four 3D MCNP input files (FeIn.mcp, NbNiAu.mcp,
  ZrAlMn.mcp, mcnp_tld.inp). Nine different foil reactions measured at 4 depths.

---

## 6. FNG-ITER Dose Rate Experiment -- DISCARDED

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
- **Evaluation:** DISCARDED. While the geometry specs are excellent (100x100 cm front,
  71.83 cm thick, cavity 119.8x150x126 mm, channel 27.4 mm ID), this benchmark
  requires a two-step R2S (Rigorous 2-Step) shutdown dose rate calculation coupling
  OpenMC neutron transport with an activation/inventory code (FISPACT). This is not
  a straightforward fixed-source transport problem -- the primary validation quantity
  is dose rate vs. cooling time, which requires depletion/activation capabilities
  beyond standard OpenMC fixed-source transport. Better tackled later after simpler
  benchmarks are established.

---

## 7. FNS Clean Benchmark: Tungsten, Vanadium, and Beryllium Assemblies -- FINALIST

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
- **Evaluation:** FINALIST. SINBAD pages provide detailed specs: W cylinder
  (629 mm diameter, 507 mm height, bricks 50.7-50.8 mm thick), V cube
  (25.4 cm side with 50.8 mm graphite reflector), source at 200 mm from assembly.
  Detector positions at 0, 76, 228, 380 mm depths. MCNP input files (mcnp-w.inp,
  mcnp-v.inp). D-T source at 350 keV deuteron beam, ~3.7E11 Bq tritium target.
  Five dosimetry reactions. Already validated with OpenMC -- highest-confidence
  starting point.

---

## 8. FNS Dogleg Duct Streaming Experiment -- FINALIST

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
- **Evaluation:** FINALIST. SINBAD page provides iron assembly dimensions
  (1700x1400x1800 mm), duct cross-section (300x300 mm), three duct leg lengths
  (1150, 600, 650 mm), NE213 detector (40 mm sphere), source yield ~4E12 n/s.
  Two complete MCNP input files (mcnp-F2.inp, mcnp-J33.inp). Four measurement
  positions with spectra and dosimetry. Excellent for weight window validation.

---

## 9. KANT Beryllium Spherical Shell Experiment -- FINALIST

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
- **Evaluation:** FINALIST. SINBAD page provides shell geometry (10 cm inner
  diameter, thicknesses 5/10/17 cm, outer diameter up to 44 cm), deuteron beam
  at 150 keV, detector at 60 degrees. MCNP model (kantmcnp.i) modifiable for all
  three thicknesses. 178-group tabulated spectra in three files. Simple spherical
  geometry with single material -- easy to build. Critical for Be neutron
  multiplication validation in breeding blanket design.

---

## 10. OKTAVIAN Sphere Experiments (Iron, Nickel, Aluminium, Silicon, Tungsten, Manganese) -- FINALIST

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
- **Evaluation:** FINALIST. SINBAD pages provide Fe sphere radius (50.32 cm, 98.69%
  Fe), Ni sphere diameter (32 cm, 99.63% Ni purity), 14.1 MeV source from 245 keV
  Cockcroft-Walton accelerator, detector at 9.5 m. MCNP models available (FE2d.i,
  NI2d.i). Already implemented in openmc_fusion_benchmarks. Simplest possible
  geometry -- ideal starting point for validating the OpenMC model-building workflow.

---

## 11. IPPE Iron Spherical Shell Transmission Experiment -- FINALIST

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
- **Evaluation:** FINALIST. SINBAD page provides all five shell geometries:
  Shell 1 (R=4.5 cm, t=2.5 cm), Shell 2 (R=12 cm, t=7.5 cm), Shell 3 (R=12 cm,
  t=10 cm), Shell 4 (R=20 cm, t=18.1 cm), Shell 5 (R=30 cm, t=28 cm). D-T source
  at 280 keV max deuteron energy, 5 mm beam spot, 2.5 ns pulse width, detector at
  6.8 m at 8 degrees. Five MCNP-4C input files (mcnp_fe1-5.inp). Systematic
  thickness variation ideal for deep-penetration transport validation.

---

## 12. FNG WCLL (Water-Cooled Lithium Lead) Mock-up Experiment -- DISCARDED

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
- **Notes:** Most recent FNG blanket experiment (2019). Tests the WCLL concept which
  is one of two EU DEMO blanket candidates. Excellent for tritium breeding validation.
- **Evaluation:** DISCARDED. The specification is a ScienceDirect journal article that
  returned HTTP 403 (paywall). No SINBAD entry with freely accessible geometry/material
  specs was found. Without access to the paper, we cannot extract the exact LiPb brick
  dimensions, EUROFER plate thicknesses, or material compositions needed to build an
  OpenMC model. Could be revisited if the paper or SINBAD entry becomes accessible.

---

## 13. EU DEMO HCPB Breeding Blanket Sector Model -- DISCARDED

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
- **Evaluation:** DISCARDED. Three disqualifying factors: (1) The specification paper
  is behind the ScienceDirect paywall (HTTP 403). (2) This is a computational-only
  benchmark with no experimental measurements for validation. (3) The complexity is
  very high -- a full tokamak sector with heterogeneous blanket modules would require
  thousands of surfaces and cells, far exceeding tractable complexity for initial
  benchmark implementation.

---

## 14. TUD Iron Slab Experiment -- FINALIST

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
- **Evaluation:** FINALIST. SINBAD page provides iron slab dimensions (100x100 cm
  front, 30 cm thick, building units 20x10x5 cm), three gap configurations (A0
  solid, A1 gap at 10 cm, A2 gap at 20 cm), source-to-slab distance 19 cm,
  slab-to-detector 300 cm, 74-degree beam angle. MCNP-4A input (MCNP.DAT, 11.8 KB).
  NE213 scintillator plus three H-proportional counters plus stilbene detector.
  17 detailed data tables in TUFE-EXP.HTM. Simple geometry, comprehensive spectral data.

---

## 15. Juelich Lithium Metal Blanket Experiment -- DISCARDED

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
- **Evaluation:** DISCARDED. The SINBAD page provides only limited geometry data
  (radial positions at 0.055-0.689 m). Critical limitation: "original Morse input
  is not available" and MCNP input details are only in a JPG figure (Figure 7) that
  cannot be parsed for exact dimensions. Without machine-readable geometry specs or
  MCNP input files, building an accurate model would require manual extraction from
  figures -- error-prone and time-consuming. The experiment dates from 1976-1984
  with less precise documentation than modern FNG/FNS experiments.

---

## 16. ORNL 14 MeV SS/Borated Polyethylene Slab Experiment -- DISCARDED

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
- **Evaluation:** DISCARDED. No dedicated SINBAD abstract page exists -- the
  specification URL just points to the SINBAD master index. The only listing found
  is "ORNL 14-MeV Neutron SS/Borated Poly Slab" under SB5_FUS with year 1979.
  No geometry dimensions, material compositions, or MCNP input files are accessible
  online. Full data requires RSICC access (rsic@ornl.gov). Insufficient freely
  available documentation to build a model without RSICC distribution.

---

## Summary Table

| # | Name | Type | Complexity | Status | Reason |
|---|------|------|------------|--------|--------|
| 1 | FNG-ITER Bulk Shield | Shield mock-up | Moderate | **FINALIST** | Full SINBAD specs, MCNP input, 14+17 measurement positions |
| 2 | FNG-ITER Streaming | Streaming/shielding | Moderate-High | **FINALIST** | 6 MCNP files, already validated with OpenMC, weight window test |
| 3 | FNG HCPB TBM Mock-up | Breeding blanket | Moderate | **FINALIST** | Exact dimensions/densities, MCNP-4C model, tritium breeding |
| 4 | FNG Copper | Material shielding | Low-Moderate | **FINALIST** | In SINBAD, already in OpenMC, simple 7-plate geometry |
| 5 | FNG Tungsten | Material shielding | Low-Moderate | **FINALIST** | DENSIMET compositions, 4 MCNP files, 9 foil reactions |
| 6 | FNG-ITER Dose Rate | Shutdown dose | Moderate | **DISCARDED** | Requires R2S coupling (not simple fixed-source) |
| 7 | FNS Clean (W/V/Be) | Material neutronics | Low | **FINALIST** | Already validated with OpenMC, MCNP inputs, comprehensive data |
| 8 | FNS Dogleg Duct | Streaming | Moderate | **FINALIST** | Full dimensions, 2 MCNP files, weight window test |
| 9 | KANT Be Shells | Neutron multiplication | Low | **FINALIST** | Spherical geometry, MCNP model, 178-group spectra |
| 10 | OKTAVIAN Spheres | Leakage spectra | Very Low | **FINALIST** | Simplest geometry, already in openmc_fusion_benchmarks |
| 11 | IPPE Iron Shells | Transmission | Very Low | **FINALIST** | 5 MCNP files, systematic thickness study |
| 12 | FNG WCLL Mock-up | Breeding blanket | Moderate | **DISCARDED** | Specs behind ScienceDirect paywall |
| 13 | EU DEMO HCPB Sector | Tokamak sector | High | **DISCARDED** | Paywall, computational-only, too complex |
| 14 | TUD Iron Slab | Material shielding | Low | **FINALIST** | MCNP input, 17 data tables, 3 gap configurations |
| 15 | Juelich Li Blanket | Breeding blanket | Moderate | **DISCARDED** | No MCNP input, geometry only in JPG figures |
| 16 | ORNL SS/BPoly Slab | Shielding | Low | **DISCARDED** | No online specs, requires RSICC access |

**Final count: 11 FINALISTS, 5 DISCARDED**

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
