# Fusion Benchmark Finalists -- Detailed Specifications for OpenMC Implementation

Selected from 16 candidates after reviewing SINBAD specification documents.
All finalists have: freely accessible geometry specs, MCNP input files,
measured experimental data for validation, and tractable complexity.

---

## 1. OKTAVIAN Sphere Experiments (#10)

**Priority:** Phase 1 (simplest geometry, already implemented)
**Complexity:** Very low -- single sphere, central point source

### Geometry
- Iron sphere: radius 50.32 cm (carbon steel, 98.69% Fe)
- Nickel sphere: diameter 32 cm (99.63% Ni purity)
- Additional materials available: Al, Si, W, Mn (six spheres total)

### Source Definition
- 14.1 MeV D-T neutrons (Cockcroft-Walton accelerator, 245 keV deuterons)
- Central point source (isotropic, angle-dependent intensity correction needed)
- ~5% lower-energy component alongside 14.1 MeV peak

### Materials
- Iron sphere: carbon steel plates, 98.69% Fe composition
- Nickel sphere: 99.63% purity nickel

### Detectors / Validation Data
- NE-213 scintillator (12.7 cm dia, 5.08 cm thick): 0.07--14 MeV neutrons
- Li-6 glass scintillator (12.7 cm dia, 2.54 cm thick): 0.01--1 MeV
- Detector distance: ~9.5 m from sphere center, 17.28-degree solid angle
- Neutron leakage spectra: 10 keV to 14 MeV via time-of-flight

### MCNP Input Files
- FE2d.i (2,618 bytes) -- recommended 2D model for iron
- NI2d.i -- recommended model for nickel
- Also: mcnp1d.inp, mcnp3d.inp (obsolete versions)

### Key References
- SINBAD Fe: https://www.oecd-nea.org/science/wprs/shielding/sinbad/oktav_fe/okfe-abs.htm
- SINBAD Ni: https://www.oecd-nea.org/science/wprs/shielding/sinbad/oktav_ni/okni-abs.htm
- Already in openmc_fusion_benchmarks repository

---

## 2. IPPE Iron Spherical Shell Transmission (#11)

**Priority:** Phase 1 (simple geometry, systematic thickness study)
**Complexity:** Very low -- spherical shells, single material, 5 variants

### Geometry (Five Shells)
| Shell | Radius (cm) | Wall Thickness (cm) |
|-------|-------------|---------------------|
| 1     | 4.5         | 2.5                 |
| 2     | 12.0        | 7.5                 |
| 3     | 12.0        | 10.0                |
| 4     | 20.0        | 18.1                |
| 5     | 30.0        | 28.0                |

### Source Definition
- 14.1 MeV D-T source (280 keV max deuteron energy)
- Ti-T target on Cu radiator (0.8 mm thick, 11 mm diameter)
- Beam spot diameter: 5 mm
- Ion pulse width: 2.5 ns, repetition period multiples of 200 ns
- Mean beam current: 1 uA (at 800 ns period)

### Materials
- Iron: atom density calculated from measured sphere weights (Table 2 in spec)

### Detectors / Validation Data
- Paraterphenyl fast scintillator: 5 cm diameter x 5 cm height
- Detection angle: 8 degrees from beam axis
- Flight path: 6.8 m
- Neutron leakage spectra: 50 keV to 15 MeV for all 5 thicknesses

### MCNP Input Files
- mcnp_fe1.inp through mcnp_fe5.inp (~6 KB each, MCNP-4C)

### Key Reference
- SINBAD: https://www.oecd-nea.org/science/wprs/shielding/sinbad/ippe-fe/ippe_fe-a.htm

---

## 3. KANT Beryllium Spherical Shells (#9)

**Priority:** Phase 1 (simple geometry, critical for breeding blanket validation)
**Complexity:** Low -- spherical shells, single material, 3 thickness variants

### Geometry
- Inner shell diameter: 10 cm (all configurations)
- Shell thicknesses: 5 cm, 10 cm, 17 cm
- Maximum outer diameter: 44 cm (17 cm shell case)

### Source Definition
- T(d,n)He-4 reaction, 150 keV deuteron beam
- Ti-T target on Cu backing
- Source strength monitored by associated alpha particles (Si surface barrier detector)

### Materials
- Beryllium metal (density ~1.85 g/cm3)

### Detectors / Validation Data
- Measurement angle: 60 degrees from deuteron beam
- NE-213 scintillator: energies above 3 MeV
- Proton recoil proportional counters: down to ~50 keV
- Time-of-flight: below 100 keV
- Bonner sphere spectrometer: source energy to thermal
- Output: 178-group neutron leakage spectra (50 keV to 15 MeV)
- Five-group partial leakage multiplications
- Combined uncertainty: 8% (1-sigma)

### MCNP Input Files
- kantmcnp.i -- model for 17 cm shell (adjust outer radius for 5/10 cm)
- fzkbe.tb1, fzkbe.tb2, fzkbe.tb3 -- 178-group tabulated spectra
- kant.xls -- measured spectra spreadsheet

### Key Reference
- SINBAD: https://www.oecd-nea.org/science//wprs/shielding/sinbad/kant/fzk-be_a.htm

---

## 4. FNS Clean Benchmarks: Tungsten, Vanadium, Beryllium (#7)

**Priority:** Phase 1 (already validated with OpenMC, comprehensive data)
**Complexity:** Low -- simple cylindrical/cubic geometry, single material per assembly

### Geometry
- **Tungsten cylinder:** diameter 629 mm, height 507 mm, brick thickness 50.7-50.8 mm
- **Vanadium cube:** 25.4 cm side, with 50.8 mm graphite reflector on 4 sides + rear
- **Source distance:** 200 mm from assembly front surface
- **Detector positions:** 0, 76, 228, 380 mm depth (along central axis)

### Source Definition
- D-T source: Ti-T target (~3.7E11 Bq activity), 350 keV deuteron beam
- Approximately isotropic point source at 14 MeV
- Angle-dependent corrections available
- MCNP source cards: si1/sp1 energy bins, sb2/vec/dir biasing cards
- Source weight factor: 1.1261

### Materials
- **Tungsten:** W-Ni-Cu alloy (exact isotopic composition in Table 1 of fnsw-exp.htm)
- **Vanadium:** >99.7% purity (composition in Table 1 of fnsv-exp.htm)
- **Graphite reflector:** composition detailed in Table 1

### Detectors / Validation Data
- NE213 scintillator (14 mm sphere in 22 mm channel): >2 MeV neutrons
- Proton recoil counter (19 mm OD, 127 mm effective length): 20 keV--1 MeV
- BF3 counter (14 mm OD, 99 mm, 96% B-10 enriched): 1--300 eV
- BC537 scintillator (40 mm sphere): gamma-ray spectra
- Dosimetry foils: Al-27(n,a), Nb-93(n,2n), In-115(n,n'), W-186(n,g), Au-197(n,g)
- TLD dosimeters: Mg2SiO4, Sr2SiO4, Ba2SiO4 for gamma heating

### MCNP Input Files
- mcnp-w.inp (12,010 bytes) -- tungsten model
- mcnp-v.inp (13,054 bytes) -- vanadium model

### Key References
- SINBAD W: https://www.oecd-nea.org/science/wprs/shielding/sinbad/fns_w/fnsw-abs.htm
- SINBAD V: https://www.oecd-nea.org/science/wprs/shielding/sinbad/fns_v/fnsv-abs.htm
- OpenMC validation: https://www.tandfonline.com/doi/full/10.1080/15361055.2024.2323747
- Already in openmc_fusion_benchmarks repository

---

## 5. FNG Copper Benchmark (#4)

**Priority:** Phase 2 (simple geometry, already benchmarked with OpenMC)
**Complexity:** Low-moderate -- 7-plate slab, pure material

### Geometry
- Assembly: 60 x 60 x 69.9 cm (seven 10-cm copper plates)
- Source distance: 5.3 cm from front surface (standard FNG geometry)

### Source Definition
- 14 MeV D-T source at FNG (Frascati Neutron Generator)
- FORTRAN source subroutine available in SINBAD

### Materials
- Pure copper (OFHC grade, exact composition in SINBAD entry)

### Detectors / Validation Data
- 7 foil reaction types: Nb-93(n,2n), Al-27(n,a), Ni-58(n,p), In-115(n,n'),
  Au-197(n,2n), W-186(n,g), Au-197(n,g)
- NE213 neutron flux spectra
- TLD dose rates
- Measured at multiple penetration depths through 7 plates

### MCNP Input Files
- Full 3D model and simplified model available in SINBAD
- Already converted to OpenMC in published studies

### Key References
- Paper: https://www.sciencedirect.com/science/article/abs/pii/S0920379618308184
  (paywall, but data in SINBAD)
- Available in openmc_fusion_benchmarks repository

---

## 6. FNG Tungsten Experiment (#5)

**Priority:** Phase 2 (important for ITER divertor design)
**Complexity:** Low-moderate -- block geometry, well-characterized alloy

### Geometry
- Block dimensions: 42-47 cm (length) x 46.85 cm (height) x 49 cm (thickness)
- DENSIMET-180 layer: 7 cm height
- Lateral access channels: 5.2 cm diameter
- Foil slots: 4.4 mm width
- Source distance: 5.3 cm from block surface

### Source Definition
- 14 MeV D-T FNG source
- Angular/energy distributions in Figures 1-2 and Tables 1-2
- FORTRAN source routine available

### Materials
- **DENSIMET-176** (~1.5 ton): 92.3% W, 2.6% Fe, 4.2% Ni
- **DENSIMET-180** (~0.25 ton): 95.0% W, 1.6% Fe, 3.4% Ni
- TLD-300 (CaF2:Tm) chips: 3.2 x 3.2 x 0.9 mm3 in 1 mm Perspex holders

### Detectors / Validation Data
- 4 depth positions from block surface
- 9 foil reactions: Au-197(n,g), Mn-55(n,g), In-115(n,n'), Ni-58(n,p),
  Fe-56(n,p), Al-27(n,a), Ni-58(n,2n), Zr-90(n,2n), Nb-93(n,2n)
- TLD-300 dosimeters: 7 chips per position, calibration range 50 mGy--4 Gy

### MCNP Input Files
- FeIn.mcp (89,670 bytes) -- 3D model for Fe/In foil reactions
- NbNiAu.mcp (91,581 bytes) -- 3D model for Nb/Ni/Au reactions
- ZrAlMn.mcp (91,649 bytes) -- 3D model for Zr/Al/Mn reactions
- mcnp_tld.inp (132,050 bytes) -- 3D model for TLD heating

### Key Reference
- SINBAD: https://www.oecd-nea.org/science/wprs/shielding/sinbad/fng_w/fngw-a.htm

---

## 7. TUD Iron Slab Experiment (#14)

**Priority:** Phase 2 (deep-penetration benchmark, well-documented)
**Complexity:** Low -- slab geometry, single material, 3 configurations

### Geometry
- Iron slab: 100 x 100 cm front face, 30 cm thick
- Building units: 20 x 10 x 5 cm iron blocks
- Three configurations:
  - A0: solid (no gap)
  - A1: 5 cm vertical gap at 10 cm from center
  - A2: 5 cm vertical gap at 20 cm from center
- Source-to-slab distance: 19 cm
- Slab-to-detector distance: 300 cm
- Source-to-detector distance: 349 cm
- Beam angle: 74 degrees (d-beam to source-slab axis)

### Source Definition
- 14 MeV D-T pulsed neutron generator
- Time distribution: proportional to exp[-(t/1.4 ns)^2]
- Angular intensity/energy distributions: Figures 1-2

### Materials
- Iron: chemical composition in TUFE-EXP.HTM tables

### Detectors / Validation Data
- NE213 scintillator: E > 1 MeV (neutrons), E > 0.2 MeV (photons)
- Three H-filled proportional detectors + stilbene scintillator: 30 keV--2.3 MeV
- Neutron and photon spectral flux measurements
- Time-of-arrival spectra with pulse-shape discrimination
- 17 comprehensive data tables in TUFE-EXP.HTM

### MCNP Input Files
- MCNP.DAT (11,833 bytes) -- 3D model for MCNP-4A

### Key Reference
- SINBAD: https://www.oecd-nea.org/science/wprs/shielding/sinbad/tud_fe/tufe-abs.htm

---

## 8. FNG HCPB Tritium Breeder Module Mock-up (#3)

**Priority:** Phase 2 (key breeding blanket benchmark)
**Complexity:** Moderate -- box geometry with layered internals

### Geometry
- **Main box (AISI-303 SS):** 31.0 x 29.0 x 30.9 cm external, 0.5 cm wall thickness
- **Rear box (AISI-316 SS):** 31.0 x 14.8 x 30.9 cm external, 0.5 cm wall thickness
- **Breeder layers:** 1.2 cm thick, separated by 1 mm SS walls
- **Source distance:** 5.3 cm from FNG target to front of assembly

### Source Definition
- 14 MeV D-T FNG source
- D-T source routines available for MCNP5/MCNPX

### Materials
- **Beryllium metal:** density 1.85 g/cm3
- **Li2CO3 powder (front section):** density 1.123 g/cm3
- **Li2CO3 powder (rear cassette):** density 0.9413 g/cm3, total mass 11690.4 +/- 0.1 g
- **Lithium isotopics:** 7.5 at% Li-6, 92.5 at% Li-7 (natural)
- **AISI-303 SS:** density 7.954 g/cm3
- **AISI-316 SS:** rear cassette structure

### Detectors / Validation Data
- Tritium production rates: Li-6(n,t) and Li-7(n,t) at 16 positions
- Neutron reaction rates: Au, Ni, Al, Nb foils at 4 positions in beryllium
- Nuclear heating: TLD-300 in breeder layers
- Four main depth positions: y = 4.2, 10.5, 16.8, 23.1 cm from block surface

### MCNP Input Files
- mcnp-hcpb.i (84,053 bytes) -- 3D model for MCNP-4C
- D-T source routines for MCNP5 and MCNPX

### Key Reference
- SINBAD: https://www.oecd-nea.org/science/wprs/shielding/sinbad/fng_hcpb/fnghcpb-a.htm

---

## 9. FNG-ITER Bulk Shield Mock-up (#1)

**Priority:** Phase 3 (multi-material ITER mock-up)
**Complexity:** Moderate -- slab geometry with multiple material layers, ~100 cells

### Geometry
- Total mock-up thickness: 94 cm
- Source-to-mockup distance: 5.3 cm
- Source aperture: 60-degree spherical cap
- Cu, SS316, and Perspex (C5O2H8) sandwich layers
- Rear section: alternating Cu/SS plates simulating magnet

### Source Definition
- 14 MeV D-T FNG source
- Source measurement accuracy: +/- 2%
- FORTRAN source subroutine available

### Materials
- Copper (first wall / magnet layers)
- Stainless steel SS316 (structural)
- Perspex / C5O2H8 (water-equivalent blanket substitute)

### Detectors / Validation Data
- **Activation foils at 14 depths (cm from front):**
  3.43, 10.32, 17.15, 23.95, 30.80, 41.85, 46.85, 53.80, 60.55, 67.40,
  74.40, 81.10, 87.75, 92.15
- **Reactions:** Nb-93(n,2n), Al-27(n,a), Ni-58(n,p), Au-197(n,g)
- **TLD-300 at 17 positions (cm from front):**
  3.36, 10.23, 17.09, 23.95, 30.81, 41.88, 46.88, 53.89, 60.72, 67.58,
  74.57, 81.20, 88.01, 95.36, 97.36, 99.76, 101.96

### MCNP Input Files
- MCNP-4A geometry input (GEOM.DAT)
- Detector positions (DETEC.DAT)

### Key Reference
- SINBAD: https://www.oecd-nea.org/science/wprs/shielding/sinbad/FNG_BLKT/FNGBKT_A.HTM

---

## 10. FNG-ITER Streaming Experiment (#2)

**Priority:** Phase 3 (weight window / variance reduction test)
**Complexity:** Moderate-high -- 3D geometry with streaming channel and cavity

### Geometry
- **Front cross-section:** 100 x 100 cm
- **Total thickness:** 94.26 cm
- **Front Cu layer:** 1 cm thick
- **Streaming channel:** 28 mm inner diameter, 39.07 cm length, 1 mm SS316 wall
- **Detector cavity:** 52 mm (z) x 148 mm (x) x 48 mm (y)
- **Rear coil block:** 47 x 47 cm, 30.8 cm depth (Cu + SS316 plates)
- **Polythene shield:** last 30 cm
- **Source distance:** 5.3 cm from mock-up (on-axis); 5.3 cm lateral shift (off-axis)

### Source Definition
- 14 MeV D-T FNG source
- Two configurations: on-axis and off-axis (channel mouth at pi/4 to beam)
- FORTRAN source subroutine: source.for (45,178 bytes)

### Materials
- Stainless steel AISI-316 (vacuum vessel, channel walls, detector box)
- Copper (front layer, rear coil block)
- Perspex (water-equivalent shielding)
- Polythene (background shielding)

### Detectors / Validation Data
- **Foil diameter:** 18 mm; thickness 1-3 mm depending on depth (0.05 mm for Au)
- **Channel positions (cm from surface):** 0.25, 12.95, 25.95, 38.65
- **Cavity positions (cm):** 39.12--43.82 (11 foils)
- **Behind-cavity positions (cm):** 46.35, 53.30, 60.05, 66.90, 73.90, 80.60, 87.25, 91.65
- **Reactions:** Nb-93(n,2n) [10.8 MeV], Al-27(n,a) [8.5 MeV],
  Ni-58(n,p) [2.9 MeV], Au-197(n,g) [thermal]
- **TLD-300:** CaF2:Tm, 0.32 x 0.32 x 0.09 cm3

### MCNP Input Files
- mcnpfoil.inp (73,096 bytes) -- 3D foil reaction rate model
- mcnp_nh.inp (86,865 bytes) -- 3D nuclear heating model
- mcnp_hss.inp (11,787 bytes) -- simplified SS heating model
- mcnp_hcu.inp (11,837 bytes) -- simplified Cu heating model
- mcnp_tld.inp (12,110 bytes) -- simplified TLD heating model
- source.for (45,178 bytes) -- FORTRAN source routine

### Key Reference
- SINBAD: https://www.oecd-nea.org/science/wprs/shielding/sinbad/fng_str/fngstr-a.htm
- OpenMC validation: https://www.tandfonline.com/doi/full/10.1080/15361055.2024.2400762

---

## 11. FNS Dogleg Duct Streaming Experiment (#8)

**Priority:** Phase 3 (streaming problem, weight window validation)
**Complexity:** Moderate -- iron slab with doubly-bent duct

### Geometry
- **Iron assembly:** 1700 x 1400 x 1800 mm
- **Duct cross-section:** 300 x 300 mm
- **Duct legs:** 1150 mm (1st horizontal), 600 mm (vertical), 650 mm (2nd horizontal)
- **Duct connections:** right angles between legs

### Source Definition
- D-T source at FNS facility, ~4.0E12 n/s yield at full beam current
- Source aligned with first horizontal duct leg

### Materials
- Iron (assembly bulk material)

### Detectors / Validation Data
- NE213 spherical scintillator (40 mm diameter): >2 MeV neutrons
- 4 measurement positions (#3, #5, #7, #9) within and beyond duct
- Dosimetry: Nb-93(n,2n), In-115(n,n'), Au-197(n,g)

### MCNP Input Files
- mcnp-F2.inp (FENDL/2 library version)
- mcnp-J33.inp (JENDL-3.3 library version)

### Key Reference
- SINBAD: https://www.oecd-nea.org/science/wprs/shielding/sinbad/fns_duct/fnsstr-a.htm

---

## Summary of Implementation Order

| Phase | Benchmark | Geometry Type | Surface Count (est.) | Key Validation |
|-------|-----------|---------------|---------------------|----------------|
| 1 | OKTAVIAN spheres | Single sphere | ~5 | Leakage spectra |
| 1 | IPPE iron shells | Spherical shells (x5) | ~10 each | Transmission spectra |
| 1 | KANT Be shells | Spherical shells (x3) | ~10 each | Multiplication spectra |
| 1 | FNS W/V/Be | Cylinder/cube | ~50 | Spectra + reactions |
| 2 | FNG copper | 7-plate slab | ~30 | 7 reaction types |
| 2 | FNG tungsten | Block assembly | ~50 | 9 reaction types |
| 2 | TUD iron slab | Slab + gaps | ~30 | Spectral flux |
| 2 | FNG HCPB TBM | Box + layers | ~100 | Tritium production |
| 3 | FNG-ITER bulk shield | Multi-layer slab | ~100 | Reactions + heating |
| 3 | FNG-ITER streaming | 3D + channel/cavity | ~200 | Streaming + reactions |
| 3 | FNS dogleg duct | Slab + bent duct | ~100 | Streaming + dosimetry |

## Discarded Candidates Summary

| # | Name | Discard Reason |
|---|------|----------------|
| 6 | FNG-ITER Dose Rate | Requires R2S activation coupling, not simple fixed-source |
| 12 | FNG WCLL Mock-up | Specs behind ScienceDirect paywall, no SINBAD entry found |
| 13 | EU DEMO HCPB Sector | Paywall + computational-only (no experiment) + too complex |
| 15 | Juelich Li Blanket | No MCNP input available, geometry only in unreadable JPG |
| 16 | ORNL SS/BPoly Slab | No online specs, requires RSICC distribution access |
