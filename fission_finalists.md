# Fission Reactor Benchmark Finalists for OpenMC Implementation

Selected from 15 candidates on 2026-03-12. Each finalist has verified publicly accessible
specifications with extractable numeric data, reference results, and tractable complexity.

---

## 1. VERA Core Physics Benchmark (Problems 1-3)

**Type:** PWR pin cell / fuel assembly / 3D assembly (Watts Bar Unit 1)
**Spec Source:** CASL-U-2012-0131-004.pdf (freely available, text-extractable)

### Fuel Pin Dimensions (Table 1 of spec)
| Parameter | Value |
|-----------|-------|
| Pellet Radius | 0.4096 cm |
| Inner Clad Radius | 0.418 cm |
| Outer Clad Radius | 0.475 cm |
| Rod Pitch | 1.26 cm |
| Rod Height | 385.1 cm |
| Fuel Stack Height | 365.76 cm |
| Plenum Height | 16.0 cm |
| End Plug Heights (x2) | 1.67 cm |
| Fuel Material | UO2 |
| Cladding Material | Zircaloy-4 |

### Guide/Instrument Tubes (Table 2)
| Parameter | Value |
|-----------|-------|
| Inner Guide Tube Radius | 0.561 cm |
| Outer Guide Tube Radius | 0.602 cm |
| Inner Instrument Tube Radius | 0.559 cm |
| Outer Instrument Tube Radius | 0.605 cm |
| Assembly Pitch | 21.50 cm |
| Inter-Assembly Half Gap | 0.04 cm |

### Material Compositions
| Material | Property | Value |
|----------|----------|-------|
| UO2 Fuel | Density | 10.257 g/cc (94.5% TD with dish/chamfer correction) |
| UO2 Fuel | Enrichments | 2.11%, 2.619%, 3.10% (three core regions) |
| Zircaloy-4 | Density | 6.56 g/cc |
| SS-304 | Density | 8.00 g/cc |
| Inconel-718 | Density | 8.19 g/cc |
| Moderator (HZP) | Density | 0.743 g/cc at 565K, 2250 psi |
| Moderator | Boron | 1300 ppm (ARO critical) |

### Isotopic Number Densities (Problem 1A, 3.1% enriched UO2, atoms/bn-cm)
| Isotope | Number Density |
|---------|---------------|
| U-234 | 6.11864E-06 |
| U-235 | 7.18132E-04 |
| U-236 | 3.29861E-06 |
| U-238 | 2.21546E-02 |
| O-16 (fuel) | 4.57642E-02 |
| H-1 (moderator) | 4.96224E-02 |
| O-16 (moderator) | 2.48112E-02 |
| B-10 (1300 ppm) | 1.07070E-05 |
| B-11 (1300 ppm) | 4.30971E-05 |

### Reference k-eff Values (CE KENO-VI, ENDF/B-VII.0)
| Problem | Conditions | k-eff |
|---------|-----------|-------|
| 1A | 565K, 0.743 g/cc mod | 1.187038 +/- 0.000054 |
| 1B | 600K, 0.661 g/cc mod | 1.182149 +/- 0.000068 |
| 1C | 900K fuel, 565K mod | 1.171722 +/- 0.000072 |
| 1D | 1200K fuel, 565K mod | 1.162603 +/- 0.000071 |
| 1E | IFBA pin, 600K | 0.771691 +/- 0.000076 |

### Geometry Layout
- 17x17 lattice: 264 fuel rods, 24 guide tubes, 1 instrument tube
- Problem 1: single pin cell (reflective BC)
- Problem 2: single 2D assembly (17 sub-cases: A-Q with varying temperatures, poisons, control rods)
- Problem 3: 3D assembly with axial detail (nozzles, spacer grids, plenum)
- Problems 4-10: multi-assembly up to full core

### Where to Find in PDF
- Geometry: Section 1 (pages 4-17), Tables 1-16
- Materials: Section 2 (page 18), Table 17 (enrichment equations)
- Operating conditions: Section 3 (page 19), Table 18
- Problem 1 specs: pages 20-23, Tables P1-1 through P1-5
- Problem 2 specs: pages 24-39, Tables P2-1 through P2-5

---

## 2. BEAVRS (Single Assembly Extraction)

**Type:** PWR fuel assembly from commercial plant (Westinghouse 4-loop, 3411 MWth)
**Spec Source:** https://github.com/mit-crpg/BEAVRS (Python API specs)

### Key Dimensions
| Parameter | Value |
|-----------|-------|
| Fuel Pellet Radius | 0.39218 cm |
| Inner Clad Radius | 0.40005 cm |
| Outer Clad Radius | 0.45720 cm |
| Rod Pitch | 1.25984 cm |
| Assembly Pitch | 21.50364 cm |
| Active Fuel Height | 365.76 cm |
| Fuel Material | UO2 |
| Cladding | Zircaloy-4 |

### Material Compositions
- Enrichments: 1.6%, 2.4%, 3.1% (three assembly types)
- UO2 density: 10.257 g/cc
- Moderator: borated water at various concentrations
- Burnable absorbers: Pyrex (borosilicate), IFBA, WABA, Gadolinia

### Reference Results
- Measured HZP critical boron concentration: ~975 ppm
- Measured boron letdown curves for 2 full cycles
- 3D in-core flux maps from 58 instrumented assemblies
- These are REAL MEASURED DATA, not computed reference solutions

### Geometry Layout
- 193 assemblies total in full core (we implement single assembly or 3x3 colorset)
- 17x17 lattice per assembly
- Baffle, barrel, neutron pad regions specified

### Implementation Note
We build from the specification data, not from the existing OpenMC model in the repo.
The full core would exceed 8GB RAM; single assembly or 3x3 colorset is the target.

---

## 3. NuScale-Like SMR (McSAFER Benchmark)

**Type:** Small modular reactor, boron-free design
**Spec Source:** MDPI paper (open access) + ECP ExaSMR GitHub

### Key Dimensions
| Parameter | Value |
|-----------|-------|
| Fuel Rod Pitch | 1.2598 cm |
| Assembly Pitch | 21.5036 cm |
| Active Fuel Length | 200 cm |
| Assembly Type | 17x17 (264 fuel rods, 24 guide tubes, 1 IT) |

### Material Compositions
- Enrichments: 1.6%, 2.4%, 3.1% UO2
- 7 assembly types with different enrichment/poison combinations
- Boron-free design (relies on burnable absorbers and control rods)

### Reference Results
- Serpent reference k-eff for multiple configurations
- Control rod worth
- Radial and axial power distributions
- OpenMC vs. Serpent: ~355 pcm difference

### Geometry Layout
- 37 fuel assemblies in full core (vs 193 for BEAVRS)
- Compact SMR geometry
- Single assembly: ~300-500 surfaces (very tractable)
- Full core: may fit in 8GB with --small-tallies

### Distinctive Feature
Boron-free design is physically distinct from VERA/BEAVRS. The shorter active fuel
length (200 cm vs 366 cm) and smaller core are characteristic of SMR designs.

---

## 4. MSRE (Molten Salt Reactor Experiment)

**Type:** Research reactor, molten salt, graphite moderated
**Spec Source:** Frontiers paper (open access) + IRPhEP MSRE-MSR-EXP-001

### Key Dimensions
| Parameter | Value |
|-----------|-------|
| Vessel Inner Diameter | 147.32 cm |
| Vessel Height | ~238.76 cm |
| Moderator | Vertical graphite stringers with fuel salt channels |

### Material Compositions
- Fuel: LiF-BeF2-ZrF4-UF4 molten salt with U-235
- Moderator: Graphite (CGB grade)
- Detailed salt composition specified in IRPhEP evaluation

### Reference Results
| Code | k-eff |
|------|-------|
| Experimental | 0.99978 +/- 0.00420 |
| OpenMC CSG | 1.00878 +/- 0.00032 |
| OpenMC CAD | 1.00872 +/- 0.00040 |
| Serpent 2.1.30 (ENDF/B-VII.1) | 1.02132 +/- 0.00003 |

### Geometry Layout
- Cylindrical vessel with graphite stringer array
- Fuel salt channels between graphite blocks
- Complex internal structure (moderator, control rods, downcomer)
- ~500-2,000 surfaces for CSG model

### Distinctive Feature
The only molten salt reactor benchmark in the list. Exercises OpenMC's capabilities
for non-standard geometries (no fuel pins). Large experimental uncertainty (420 pcm)
but good code-to-code agreement.

---

## 5. HTR-10 Pebble Bed Reactor

**Type:** High-temperature gas-cooled reactor, pebble bed, TRISO fuel
**Spec Source:** IAEA HTR-10 benchmark reports + numerous published papers

### Key Dimensions
| Parameter | Value |
|-----------|-------|
| Core Diameter | 180 cm |
| Critical Height | ~123-126 cm (depending on configuration) |
| Pebble Outer Radius | 3.0 cm |
| Fuel Zone Radius | 2.5 cm |
| Graphite Shell Thickness | 0.5 cm |

### TRISO Particle Dimensions
| Layer | Outer Radius |
|-------|-------------|
| UO2 Kernel | 0.025 cm (250 um) |
| Buffer (PyC) | 0.034 cm |
| Inner PyC | 0.038 cm |
| SiC | 0.0415 cm |
| Outer PyC | 0.0455 cm |

### Material Compositions
| Material | Property | Value |
|----------|----------|-------|
| UO2 Kernel | Enrichment | 17% U-235 |
| UO2 Kernel | Density | ~10.4 g/cc |
| Graphite Matrix | Density | ~1.73 g/cc |
| Graphite Shell | Density | ~1.73 g/cc |
| Graphite Reflector | Density | ~1.76 g/cc |
| Coolant | Material | Helium |
| TRISO per Pebble | Count | 8,335 |

### Reference Results
- Critical height measurements for multiple configurations
- VSOP, MCNP, Serpent, and OpenMC reference k-eff values
- Multiple packing fraction and homogenization studies

### Geometry Layout
- Cylindrical core with random pebble packing
- Fuel pebbles + moderator-only pebbles
- Graphite reflector surrounding core
- Uses OpenMC's `openmc.model.TRISO` for particle packing

### Distinctive Feature
Double heterogeneity (TRISO in pebbles, pebbles in core) is the signature modeling
challenge. Exercises OpenMC's unique random packing capabilities not found in
deterministic codes. Non-LWR advanced reactor type.

---

## 6. KRITZ-2 Critical Experiments

**Type:** Critical experiment, UO2 and MOX lattices at varying temperatures
**Spec Source:** OECD/NEA NEA/NSC/DOC(2005)24 (freely available PDF, text-extractable)

### Configuration Specifications (Table 2.1)
| Parameter | KRITZ-2:1 (UO2) | KRITZ-2:13 (UO2) | KRITZ-2:19 (MOX) |
|-----------|-----------------|-------------------|-------------------|
| Fuel | UO2, 1.86 wt% 235U | UO2, 1.86 wt% 235U | MOX, 1.5 wt% PuO2 |
| Pu Composition | -- | -- | 91.41 at% 239Pu |
| Fuel Diameter | 10.58 mm | 10.58 mm | 9.45 mm |
| Clad Outer Diameter | 12.25 mm | 12.25 mm | 10.79 mm |
| Rod Pitch | 14.85 mm | 16.35 mm | 18.00 mm |
| Number of Rods | 44 x 44 | 40 x 40 | 25 x 24 |
| Temp (cold) | 19.7 C | 22.1 C | 21.1 C |
| Temp (hot) | 248.5 C | 280.1 C | 235.9 C |
| Boron (cold, ppm) | 652.8 | 451.9 | 665.6 |
| Boron (hot, ppm) | 1055.2 | 961.7 | 1000.1 |
| Water Height (cold, mm) | -- | -- | -- |
| Water Height (hot, mm) | -- | -- | -- |

### Reference k-eff Values (Average of 13 solutions)
| Configuration | k-eff | Std Dev |
|--------------|-------|---------|
| KRITZ-2:1, 19.7 C | 0.99463 | 0.00330 |
| KRITZ-2:1, 248.5 C | 0.99290 | 0.00346 |
| KRITZ-2:13, 22.1 C | 0.99736 | 0.00204 |
| KRITZ-2:13, 280.1 C | 0.99605 | 0.00280 |
| KRITZ-2:19, 21.1 C | 0.99908 | 0.00179 |
| KRITZ-2:19, 235.9 C | 0.99822 | 0.00177 |

### Geometry Layout
- Rectangular lattice of fuel rods in light water
- Simple geometry (uniform pin pitch, no guide tubes or control rods)
- Full specifications in Appendix A of the OECD/NEA report

### Distinctive Feature
Temperature-dependent criticality: same core measured at room temperature and ~245 C.
Provides unique validation of Doppler reactivity and moderator density effects.
MOX core (KRITZ-2:19) provides weapons-grade Pu validation. No existing OpenMC model.

### Where to Find in PDF
- Table 2.1: Configuration summary (page 13)
- Appendix A: Full detailed specifications
- Tables 4.1-4.6: k-eff results by configuration
- Tables 4.9-4.13: Pin power distributions
- Tables 4.14-4.16: Pin cell k-inf values

---

## 7. VVER-1000 Mock-up Critical Experiments

**Type:** Hexagonal lattice PWR mock-up (LR-0 research reactor, Czech Republic)
**Spec Source:** Open access paper (Nuclear Engineering and Technology, 2021)

### Key Dimensions
| Parameter | Value |
|-----------|-------|
| Assembly Pitch | 23.6 cm |
| Fuel Pin Pitch | 12.75 mm (1.275 cm) |
| Fuel Pin Length | ~1.35 m |
| Fissile Column Length | 1.25 m |
| Number of Assemblies | 32 |
| Fuel Pins per Assembly | 312 (282 in central assembly #27) |
| Central Channel Diameter | 6.8 cm (dry, cladded with 2.5 mm Zr alloy) |

### Material Compositions
- Fuel enrichments: 2.0%, 3.0%, 3.3% 235U
- Cladding: Zirconium alloy
- Absorber clusters: 18 per assembly
- Moderator: Light water with boric acid

### Reference Results
- Six experimental critical configurations (keff = 1.0)
- MCS/MCNP6 with ENDF/B-VII.1: keff overprediction +137 to +532 pcm
- MCS/MCNP6 with ENDF/B-VIII.0: reduced overprediction (max diff 27 +/- 25 pcm vs VII.1)
- Four pin-by-pin power maps (96-260 pins per map)
- MCS and MCNP6 input files available as supplementary material

### Geometry Layout
- Hexagonal fuel assembly lattice (distinctive from all other candidates)
- Baffle, barrel, displacer, pressure vessel, biological shielding simulators
- Six critical configurations with varying water level and boric acid concentration

### Distinctive Feature
The only hexagonal-lattice benchmark in the finalist list. Exercises OpenMC's
`openmc.HexLattice` capabilities. VVER geometry is fundamentally different from
Western square-lattice PWR designs. No existing OpenMC model.

### Where to Find in Paper
- Section 2: Benchmark description with key dimensions
- Table 1: Critical configuration names
- Table 2: Critical configurations (moderator level, boric acid)
- Tables 3-4: keff results by library
- Tables 5-8: Sensitivity analysis
- Table 9: (C-E) keff discrepancies and uncertainties
- Supplementary material: MCS and MCNP6 input files

---

## 8. SNAP-10A/2 Space Reactor

**Type:** Compact space reactor, UZrH fuel, critical experiments
**Spec Source:** ORNL/TM-2005/54 (freely available PDF, text-extractable)

### Key Dimensions
| Parameter | Value |
|-----------|-------|
| Fuel Element Diameter | 1.250 in (3.175 cm) |
| Fuel Element Length | 12.450 in (31.623 cm) |
| Triangular Pitch | 1.260 in (3.2004 cm) |
| Number of Fuel Elements | 37 (full core) |
| Core Vessel Inside Diameter | 8.900 in (22.606 cm) |

### Material Compositions
| Material | Property | Value |
|----------|----------|-------|
| Fuel | Composition | 10 wt% U, 90 wt% Zr (hydride) |
| Fuel | U Enrichment | >= 93 wt% 235U |
| Fuel | Mass Density | ~6.06 g/cc |
| Fuel | H Atom Density | ~6.5E22 atoms/cm3 |
| Fuel | 235U per Element | ~128.5 g |
| Full Core | Total 235U | ~4.75 kg |
| Fuel | Burnable Poison | Trace Sm2O3 (2.7-4.0 g total per core) |

### Reference Results
- 73 benchmark configurations (56 critical, 17 subcritical)
- Various water reflection and immersion conditions
- MCNP5 reference calculations for all configurations
- Detailed reactivity measurements

### Geometry Layout
- Cylindrical stainless steel vessel
- 37 fuel elements on triangular pitch
- Up to 6 beryllium internal reflector inserts
- Various external reflector/water configurations
- ~100-300 surfaces (simplest geometry of all finalists)

### Distinctive Feature
Unique fuel type (uranium-zirconium hydride), compact space reactor geometry.
Hydrogen in fuel acts as both moderator and fuel matrix. Very different neutron
spectrum from water-moderated systems. Simple geometry allows quick implementation.

### Where to Find in PDF
- Section 3.1.1: Reactor core vessel (page 9)
- Section 3.1.2: Reflector tanks (page 11)
- Section 3.1.3: SCA-4 fuel elements (page 14)
- Table 1: Summary of 56 critical benchmarks (page 17)
- Table 2: Summary of 17 subcritical benchmarks
- Section 4: Detailed MCNP5 model descriptions

---

## 9. IAEA 10 MW MTR Research Reactor

**Type:** Pool-type materials testing reactor, plate-type fuel
**Spec Source:** IAEA-TECDOC-643 Volume 1 (freely available PDF)

### Key Specifications
| Parameter | HEU Core | LEU Core |
|-----------|----------|----------|
| Enrichment | 93% 235U | 20% 235U |
| Plates per Element | 23 | 20 (or other configs) |
| 235U per Element | 280 g | 390 g |
| Fuel Meat | UAl alloy | U3Si2-Al |
| Core Layout | 5x6 elements | 5x6 elements |
| Fuel Elements | 23 standard | 23 standard |
| Control Elements | 5 | 5 |
| Reactor Power | 10 MW | 10 MW |
| Coolant | Light water | Light water |
| Reflector | Beryllium + graphite | Beryllium + graphite |

### Reference Results
- k-eff for fresh and burned cores
- Power peaking factors (radial x local)
- Control rod worths
- Isothermal reactivity coefficients
- Transient analysis results (LOCA, rod withdrawal)
- Multiple independent calculations from international participants

### Geometry Layout
- 5x6 array of fuel elements in water pool
- Plate-type fuel (rectangular geometry, not cylindrical pins)
- Beryllium and graphite reflector regions
- ~500-1,500 surfaces

### Distinctive Feature
Plate-type fuel geometry is unique among all finalists. Exercises OpenMC's
rectangular prism modeling (vs cylindrical pins in all other benchmarks).
The IAEA MTR benchmark is arguably the most widely reproduced research reactor
benchmark in the world, with dozens of published solutions across many codes.
HEU-to-LEU conversion comparison adds validation value.

### Where to Find in PDF
- Chapter 2: Summary of benchmark analysis results
- Chapter 7: Detailed benchmark results with tables
- Appendix G: Full benchmark specifications (G-1 through G-6)
- Appendix A: Neutronics analysis details
- Table 7.1: Kinetic parameters
- Table 7.2: Reactivity coefficients
- Tables 7.5-7.8: Power peaking and control rod worths

---

## 10. KRITZ-2:19 MOX Pin Cell (Bonus -- Extracted from KRITZ-2)

This is not a separate benchmark but a note that the KRITZ-2:19 MOX configuration
can be implemented as a standalone pin cell exercise in addition to the full core
model. The MOX fuel with weapons-grade plutonium (91.41 at% 239Pu) provides
a validation case for plutonium-bearing fuel that no other finalist covers.

---

## Implementation Priority Order

Recommended implementation order based on specification quality, diversity of
reactor types, and progressive complexity:

1. **VERA Problem 1** (pin cell) -- simplest possible, validates basic setup
2. **VERA Problem 2** (2D assembly) -- moderate complexity, extensive reference data
3. **SNAP-10A/2** (space reactor) -- simple geometry, unique fuel type
4. **KRITZ-2** (critical experiments) -- temperature validation, MOX fuel
5. **BEAVRS** (single assembly) -- real measured data
6. **NuScale SMR** (single assembly) -- boron-free SMR design
7. **VVER-1000** (hex assembly) -- hexagonal lattice validation
8. **IAEA MTR** (plate fuel) -- plate-type geometry
9. **HTR-10** (pebble bed) -- TRISO double heterogeneity
10. **MSRE** (molten salt) -- unique non-LWR reactor type

### Diversity Coverage
- **Square lattice PWR:** VERA, BEAVRS, NuScale (3 benchmarks)
- **Hexagonal lattice:** VVER-1000 (1 benchmark)
- **Plate-type fuel:** IAEA MTR (1 benchmark)
- **Pebble bed HTGR:** HTR-10 (1 benchmark)
- **Molten salt:** MSRE (1 benchmark)
- **Space reactor / UZrH:** SNAP-10A/2 (1 benchmark)
- **Temperature-dependent / MOX:** KRITZ-2 (1 benchmark)

### Computational Requirements
- Pin cell (VERA P1, KRITZ-2 cell): < 1 GB RAM, < 1 min
- Single assembly (VERA P2, BEAVRS, NuScale): 1-2 GB RAM, 5-30 min
- Small core (SNAP-10A, KRITZ-2 core, VVER-1000): 2-4 GB RAM, 30-120 min
- Full MTR core (IAEA MTR): 2-4 GB RAM, 30-60 min
- Pebble bed (HTR-10): 2-6 GB RAM, 60-180 min (packing-dependent)
- Molten salt (MSRE): 2-4 GB RAM, 30-120 min
