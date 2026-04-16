# FHWA/NRC Dry Inclined Cable Tests

Primary role: damping and Scruton-number design reference for dry inclined
cable vibration.

## Source

| Ref | Citation | DOI/URL | Access status |
|---|---|---|---|
| S1 | S. Kumarasena, N.P. Jones, P. Irwin, P. Taylor, "Wind-Induced Vibration of Stay Cables", FHWA-HRT-05-083, Federal Highway Administration, 2007 | https://www.fhwa.dot.gov/publications/research/infrastructure/bridge/05083/chap3.cfm | Public FHWA report |
| S2 | FHWA-HRT-05-083 Appendix C, "Wind-Induced Cable Vibrations" | https://www.fhwa.dot.gov/publications/research/infrastructure/bridge/05083/appendc.cfm | Public FHWA report |

## Test Rig and Model

| Quantity | Value | Unit | Source | Status |
|---|---:|---|---|---|
| Cable model length | 6.7 | m | S1 | published |
| Directly exposed length | at least 5.9 | m | S1 | published |
| Diameter | 160 | mm | S1 | published |
| Effective mass per unit length | 60.8 | kg/m | S1 | published |
| Low damping range | 0.03-0.09 | percent critical | S1 | published |
| Intermediate damping range | 0.05-0.10 | percent critical | S1 | published |
| High damping range | 0.15-0.25 | percent critical | S1 | published |
| Very high damping range | 0.30-1.00 | percent critical | S1 | published |
| Surface conditions | smooth and rough PE tube | - | S1 | published |

## Test Matrix Excerpt

| Setup | Cable inclination theta | Yaw beta | Notes |
|---|---:|---:|---|
| 1B | 45 deg | 0 deg | S1 Table 1 |
| 1C | 30 deg | 35.3 deg | S1 Table 1 |
| 2A | 60 deg | 0 deg | S1 Table 1 |
| 2C | 45 deg | 45 deg | S1 Table 1 |
| 3A | 35 deg | 0 deg | S1 Table 1 |
| 3B | 20 deg | 29.4 deg | S1 Table 1 |

## Observed Response

| Quantity | Value | Unit | Source | Status |
|---|---:|---|---|---|
| Limited-amplitude vibration | up to 1D double amplitude | - | S1 | published |
| Significant vibration threshold | mainly zeta < 0.001 | - | S1 | published |
| Suppression threshold | zeta > 0.003 gives no significant vibration above 10 mm | - | S1 | published |

## Representative Input Values

| Quantity | Value | Unit | Source | Status |
|---|---:|---|---|---|
| Cable length | 6.7 | m | S1 | published |
| Cable inclination | 45 | deg | S1 setup 1B | published |
| Diameter | 0.160 | m | S1 | published |
| Mass per unit length | 60.8 | kg/m | S1 | published |
| Representative frequency | 1.0 | Hz | local model | assumed for initial tension estimate |
| Mean tension | 1.092e4 | N | local model | inferred from `T=(2 L f)^2 m` |
| EA | 4.021e9 | N | local model | inferred from circular area and E=200 GPa |

## Files

| File | Purpose |
|---|---|
| `settings.json` | Standalone GUI/C++ input for a representative dry-inclined cable |
| `cables/FHWA_Dry_1B.json` | Per-cable detail file based on setup 1B geometry |

## Modeling Notes

- This is an aerodynamic stability reference, not an RWIV case.
- `settings.json` distributes the mean tension into top and bottom endpoint
  tensions so the gravity-loaded static initial condition is internally
  consistent.
- The current input is structural-only. Aerodynamic instability must be modeled
  separately if the amplitude-growth behavior is to be reproduced.
