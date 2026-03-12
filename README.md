# Conjunction Assessment SDN Plugin

Propagator-agnostic conjunction assessment engine for the [Space Data Network](https://github.com/the-lobsternaut/space-data-network). Finds close approaches between space objects, computes collision probability using five methods, and outputs CDM FlatBuffers.

**Any propagator. Any Pc method. 47K catalog pairs screened with KD-tree + pthreads.**

## Architecture

```
┌──────────────────────────────────────────────────────────┐
│                  EphemerisSource (interface)              │
│                                                          │
│  SGP4EphemerisSource   ← TLE / GP elements               │
│  GPEphemerisSource     ← CelesTrak OMM JSON              │
│  OEMEphemerisSource    ← Tabulated ephemeris (Hermite)    │
│  CallbackEphemerisSource ← Your function pointer          │
└──────────────┬───────────────────────────────────────────┘
               │ state_at(jd)
               ▼
┌──────────────────────────────────────────────────────────┐
│                 ConjunctionEngine                        │
│                                                          │
│  1. TCA Finding     (golden-section refinement)          │
│  2. B-plane Build   (Montenbruck & Gill §6.5)           │
│  3. Mahalanobis     (2D encounter + 3D inertial)         │
│  4. Pc Computation  (pluggable — see below)              │
│  5. CDM Output      ($CDM FlatBuffer)                    │
└──────────────┬───────────────────────────────────────────┘
               │
               ▼
┌──────────────────────────────────────────────────────────┐
│                  PcMethod (interface)                     │
│                                                          │
│  AlfanoMaxPc    — Upper bound        (AAS 03-548)        │
│  Foster2D       — Short-encounter    (NASA JSC, 1992)    │
│  Patera2001     — 2D integration     (JGCD 24(4))       │
│  Chan2008       — Polar integration  (Aerospace Press)   │
│  Alfriend2D     — Gauss-Legendre     (Space Debris 1(1)) │
└──────────────────────────────────────────────────────────┘
```

## Collision Probability Methods

| Method | Type | Best For | Reference |
|---|---|---|---|
| `alfano` | Analytical bound | Quick screening, upper bound | Alfano, AAS 03-548 |
| `foster` | Short-encounter | CDM standard, most conjunctions | Foster-Estes, NASA/JSC-25898, 1992 |
| `patera` | 2D numerical | Elliptical covariance | Patera, JGCD 24(4), 2001 |
| `chan` | Polar integration | NASA CARA operations | Chan, Aerospace Press, 2008 |
| `alfriend` | Gauss-Legendre | Reference / validation | Alfriend & Akella, Space Debris 1(1), 1999 |

All methods return `PcResult` with probability, Mahalanobis distance, and Alfano upper bound.

## Mahalanobis Distance

Two variants computed for every conjunction event:

| Variant | Formula | Meaning |
|---|---|---|
| **2D (B-plane)** | √(Δr^T × C₂ₓ₂⁻¹ × Δr) | Miss in σ-units in the encounter plane |
| **3D (inertial)** | √(Δr^T × C₃ₓ₃⁻¹ × Δr) | Miss in σ-units in full 3D position space |

## Screening Engine

High-performance catalog screening for close approaches:

1. **Perigee/apogee prefilter** — eliminate altitude-incompatible pairs (O(n))
2. **KD-tree spatial indexing** — O(n log n) candidate identification
3. **Dynamic windowing** — adaptive step size based on closing rate
4. **Golden-section TCA refinement** — sub-millisecond TCA accuracy
5. **Threaded** — pthreads for native, WASM pthreads for browser

### Performance

| Catalog Size | Screening Time | Pairs Checked |
|---|---|---|
| 100 objects | ~50 ms | 4,950 |
| 1,000 objects | ~2 sec | 499,500 |
| 47,000 (full catalog) | ~minutes | 1.1B (with prefilter) |

## Usage

### Propagator-Agnostic Assessment

```cpp
#include "conjunction/conjunction_engine.h"
#include "conjunction/ephemeris_source.h"

// Any propagator → EphemerisSource
SGP4EphemerisSource primary(tle);           // from TLE
OEMEphemerisSource secondary(oem_points);   // from HPOP output
CallbackEphemerisSource custom(my_func);    // from anything

// Pick your Pc method
ConjunctionEngine engine;
engine.set_pc_method("foster");  // or "alfano", "patera", "chan", "alfriend"

// Full assessment
auto event = engine.assess(primary, secondary, start_jd, 7.0);

printf("TCA: %s\n", event.tca_iso.c_str());
printf("Miss: %.3f km\n", event.miss_distance_km);
printf("Pc (Foster): %.6e\n", event.pc.probability);
printf("Mahalanobis 2D: %.2f σ\n", event.mahalanobis_2d);
printf("Mahalanobis 3D: %.2f σ\n", event.mahalanobis_3d);
```

### Pc-Only (No Propagation)

```cpp
// Pre-computed states at TCA (from any source)
StateVector s1 = { jd, x1, y1, z1, vx1, vy1, vz1 };
StateVector s2 = { jd, x2, y2, z2, vx2, vy2, vz2 };

auto cov1 = Covariance3x3::from_rtn_diagonal(0.1, 0.3, 0.1);  // km
auto cov2 = Covariance3x3::from_rtn_diagonal(0.05, 0.15, 0.05);

auto event = engine.compute_pc(s1, s2, cov1, cov2, 0.01);  // 10m combined radius
```

### Switching Pc Methods at Runtime

```cpp
engine.set_pc_method("alfano");
auto r1 = engine.compute_pc(s1, s2, c1, c2);

engine.set_pc_method("chan");
auto r2 = engine.compute_pc(s1, s2, c1, c2);

// Compare: r1.pc.probability vs r2.pc.probability
// Alfano is always >= all other methods (it's the maximum probability)
```

### Catalog Screening

```cpp
#include "conjunction/screening.h"

ScreeningConfig config;
config.start_jd = epoch_jd;
config.duration_days = 7.0;
config.threshold_km = 5.0;
config.num_threads = 8;

ConjunctionScreener screener(config);
auto results = screener.screen(gp_elements);  // vector<GPElement>

for (auto& event : results) {
    printf("%s vs %s: %.3f km, Pc=%.2e\n",
        event.obj1_name.c_str(), event.obj2_name.c_str(),
        event.min_range_km, event.max_probability);
}
```

### GP/OMM JSON Ingestion

```cpp
#include "conjunction/gp_json.h"

// Parse CelesTrak GP JSON array
auto elements = parse_gp_json_array(json_string);

// Direct to SGP4
auto tle = gp_to_tle(gp_element);
auto state = propagate_sgp4(tle, target_jd);
```

### CDM FlatBuffer Output

```cpp
#include "conjunction/cdm_output.h"

// Assessment result → $CDM FlatBuffer
auto cdm_buffer = build_cdm_flatbuffer(event);
// cdm_buffer has $CDM file identifier
```

## Data Flow

```
Input Sources                    Engine                        Output
─────────────                    ──────                        ──────

TLE text        ─┐
GP/OMM JSON     ─┤              ┌──────────────┐
OEM ephemeris   ─┼─→ EphemerisSource ─→ │ TCA Finding    │
HPOP output     ─┤              │ B-plane Build  │ ─→ CDM ($CDM FlatBuffer)
External API    ─┘              │ Mahalanobis    │ ─→ ConjunctionEvent2 struct
                                │ Pc Computation │ ─→ JSON (via WASM embind)
                 Covariance ──→ │   (pluggable)  │
                                └──────────────┘
```

## Building

```bash
# Native build
cd src/cpp && mkdir -p build && cd build
cmake ..
make -j4

# Run tests
./test_pc_methods      # 30 tests — Pc method cross-validation
./test_engine          # 34 tests — propagator-agnostic engine
./test_gp_json         # GP/OMM JSON parsing
./test_cdm_output      # CDM FlatBuffer generation
./test_screening       # Catalog screening
```

## File Structure

```
src/cpp/
├── include/conjunction/
│   ├── conjunction_engine.h     # v2 engine (propagator-agnostic)
│   ├── ephemeris_source.h       # Propagator interface
│   ├── pc_method.h              # Pc method interface + 5 implementations
│   ├── conjunction_assessment.h # v1 engine (SGP4-coupled, backward compat)
│   ├── screening.h              # KD-tree + threaded catalog screener
│   ├── kdtree.h                 # 3D KD-tree for spatial queries
│   ├── sgp4_propagator.h        # SGP4 wrapper (dnwrnr/sgp4)
│   ├── gp_json.h                # GP/OMM JSON parser
│   └── cdm_generated.h          # CDM FlatBuffer schema
├── src/
│   ├── conjunction_engine.cpp   # TCA + B-plane + Mahalanobis
│   ├── pc_method.cpp            # 5 Pc methods
│   ├── ephemeris_source.cpp     # OEM Hermite interpolation
│   ├── conjunction_assessment.cpp # v1 (backward compat)
│   ├── screening.cpp            # Threaded screening
│   ├── kdtree.cpp               # KD-tree build + range query
│   ├── gp_json.cpp              # CelesTrak JSON parsing
│   └── cdm_output.cpp           # CDM FlatBuffer generation
└── tests/
    ├── test_pc_methods.cpp      # Pc cross-validation, Mahalanobis
    ├── test_engine.cpp          # Mixed sources, v1/v2 cross-check
    ├── test_gp_json.cpp         # GP/OMM parsing
    ├── test_cdm_output.cpp      # CDM output
    ├── test_screening.cpp       # Screening engine
    └── test_socrates.cpp        # SOCRATES-style screening
```

## References

- **Alfano, S.** "Relating Position Uncertainty to Maximum Conjunction Probability." AAS 03-548, 2003.
- **Foster, J.L. & Estes, H.S.** "A Parametric Analysis of Orbital Debris Collision Probability and Maneuver Rate for Space Vehicles." NASA/JSC-25898, 1992.
- **Patera, R.P.** "General Method for Calculating Satellite Collision Probability." JGCD 24(4), 2001.
- **Chan, F.K.** *Spacecraft Collision Probability.* El Segundo: Aerospace Press, 2008.
- **Alfriend, K.T. & Akella, M.R.** "Probability of Collision Error Analysis." Space Debris 1(1), 1999.
- **Montenbruck, O. & Gill, E.** *Satellite Orbits.* Springer, 2000. §6.5 (B-plane).
- **Vallado, D.A.** *Fundamentals of Astrodynamics and Applications.* 4th ed., 2013.

## License

Apache-2.0
