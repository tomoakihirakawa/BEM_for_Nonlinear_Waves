# Dongting Lake Bridge RWIV Field Measurements

Primary role: full-scale rain-wind-induced vibration benchmark.

## Source

| Ref | Citation | DOI/URL | Access status |
|---|---|---|---|
| S1 | Y.Q. Ni, X.Y. Wang, Z.Q. Chen, J.M. Ko, "Field observations of rain-wind-induced cable vibration in cable-stayed Dongting Lake Bridge", Journal of Wind Engineering and Industrial Aerodynamics, 95(5), 303-328, 2007 | https://doi.org/10.1016/j.jweia.2006.07.001 | Abstract and metadata public; full text may need library access |

## Bridge and Monitoring System

| Quantity | Value | Unit | Source | Status |
|---|---:|---|---|---|
| Total bridge length | 880 | m | S1 abstract/preview | published |
| Main spans | 310, 310 | m | S1 preview | published |
| Side spans | 130, 130 | m | S1 preview | published |
| Number of stay cables | 222 | - | S1 preview | published |
| Cable length range | 28-201 | m | S1 preview | published |
| Cable diameter range | 99-159 | mm | S1 preview | published |
| Field measurement duration | 45 | day | S1 abstract | published |
| Accelerometers | 15 | - | S1 preview | published |
| Displacement transducers | 1 | - | S1 preview | published |
| 3D anemometers | 2 | - | S1 preview | published |
| Rain gauges | 1 | - | S1 preview | published |

## Selected Representative Cable

The solver input in this directory uses a representative cable based on the
instrumented stay described in S1.

| Quantity | Value | Unit | Source | Status |
|---|---:|---|---|---|
| Cable name in this dataset | A12_representative | - | local | assigned |
| Cable length | 122 | m | S1 preview | published |
| Cable inclination | 30 | deg | local model | assumed, to be replaced by exact table value |
| Cable diameter | 140 | mm | S1 bridge range | assumed mid-range |
| Mass per unit length | 80 | kg/m | local model | assumed typical stay value |
| Fundamental frequency | 0.65 | Hz | local model | assumed for initial tension estimate |
| Mean tension | 2.012e6 | N | local model | inferred from `T=(2 L f1)^2 m` |
| EA | 3.079e9 | N | local model | inferred from circular area and E=200 GPa |

## Observed RWIV Response

| Quantity | Value | Unit | Source | Status |
|---|---:|---|---|---|
| Rain-wind events analyzed | 3 | - | S1 abstract | published |
| Maximum acceleration | 10 | g | S1 abstract | published |
| Maximum displacement | 0.7 | m peak-to-peak | S1 abstract | published |
| Deck-level mean wind speed | 6-14 | m/s | S1 abstract | published |
| Relative yaw angle | 10-50 | deg | S1 abstract | published |
| Rainfall | < 8 | mm/h | S1 abstract | published |
| Overall dominant mode | 3 | - | S1 abstract | published |

## Files

| File | Purpose |
|---|---|
| `settings.json` | Standalone GUI/C++ input for the representative cable |
| `cables/A12_representative.json` | Per-cable detail file with the same representative cable data |

## Modeling Notes

- The present input is a consistent simulation seed, not a calibrated
  reproduction of S1.
- `settings.json` distributes the mean tension into top and bottom endpoint
  tensions so the gravity-loaded static initial condition is internally
  consistent.
- Exact cable diameter, mass, tension, inclination, and modal frequencies
  should be replaced after full-table extraction.
- The uniform wind in `settings.json` uses 10 m/s, the midpoint of the
  published 6-14 m/s event range.
