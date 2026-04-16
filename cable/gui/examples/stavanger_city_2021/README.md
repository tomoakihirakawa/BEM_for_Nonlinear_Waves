# Stavanger City Bridge Wet/Dry Stay-Cable Observations

Primary role: open-access full-scale wet/dry comparison case.

## Source

| Ref | Citation | DOI/URL | Access status |
|---|---|---|---|
| S1 | N. Daniotti, J.B. Jakobsen, J. Snæbjörnsson, E. Cheynet, J. Wang, "Observations of bridge stay cable vibrations in dry and wet conditions: A case study", Journal of Sound and Vibration, 503, 116106, 2021 | https://doi.org/10.1016/j.jsv.2021.116106 | Open access |

## Bridge and Monitoring System

| Quantity | Value | Unit | Source | Status |
|---|---:|---|---|---|
| Bridge total length | 1067 | m | S1 | published |
| Cable-stayed suspended span | 185 | m | S1 | published |
| Tower height | 70 | m | S1 | published |
| Deck width | 15.5 | m | S1 | published |
| Deck depth | 2.4 | m | S1 | published |
| Stay cable type | locked-coil wire | - | S1 | published |
| Stay cable length range | 61-141 | m | S1 | published |
| Stay cable diameter | 79 | mm | S1 | published |
| Monitoring start | June 2019 | - | S1 | published |
| Sensors | accelerometers, anemometer, field camera | - | S1 | published |

## Instrumented Stays

| Stay | Count | Length | Diameter | Inclination | Mass | First frequency | Source | Status |
|---|---:|---:|---:|---:|---:|---:|---|---|
| C1E/C1W | 4 | 98.3 m | 79 mm | 29.8 deg | 40 kg/m | 1.05 Hz | S1 Table 1 | published |
| C2E/C2W | 4 | 95.6 m | 79 mm | 30.7 deg | 40 kg/m | 1.07-1.08 Hz | S1 Table 1 | published |

The solver input uses C1E as the representative cable.

## Observed RWIV Response

| Quantity | Value | Unit | Source | Status |
|---|---:|---|---|---|
| Critical mean wind speed | 8-12 | m/s | S1 abstract | published |
| Large-response mode | 3 | - | S1 abstract | published |
| Large-response amplitude | about 2D | peak-to-peak | S1 abstract | published |
| Reduced velocity during large response | about 35 | - | S1 abstract | published |
| Rain condition | light rainfall promotes onset | - | S1 abstract | published |

## Representative Input Values

| Quantity | Value | Unit | Source | Status |
|---|---:|---|---|---|
| Cable | C1E | - | local | assigned |
| Cable length | 98.3 | m | S1 Table 1 | published |
| Diameter | 0.079 | m | S1 Table 1 | published |
| Inclination | 29.8 | deg | S1 Table 1 | published |
| Mass per unit length | 40 | kg/m | S1 Table 1 | published |
| First frequency | 1.05 | Hz | S1 Table 1 | published |
| Mean tension | 1.705e6 | N | local | inferred from `T=(2 L f1)^2 m` |
| EA | 9.803e8 | N | local | inferred from circular area and E=200 GPa |

## Files

| File | Purpose |
|---|---|
| `settings.json` | Standalone GUI/C++ input using C1E and 10 m/s uniform wind |
| `cables/C1E.json` | Per-cable detail file using published C1E properties |

## Modeling Notes

- Rigid cross-ties between individual cables are not included yet.
- `settings.json` distributes the mean tension into top and bottom endpoint
  tensions so the gravity-loaded static initial condition is internally
  consistent.
- The cable is modeled as a single equivalent stay. The grouped-cable geometry
  in S1 should be added later if cross-tie dynamics become part of the study.

## Verification (2026-04-16)

Recomputed with:

```bash
./build_solver/cable_solver gui/examples/stavanger_city_2021/settings.json /tmp/cable_verify_stavanger_city_2021/
```

| Quantity | Published / input basis | Solver / derived value | 誤差 |
|---|---:|---:|---:|
| Cable length | 98.3 m | 98.3 m | exact |
| Line density | 40 kg/m | 40 kg/m | exact |
| Diameter | 79 mm | 79 mm | exact |
| First natural frequency f₁ | 1.05 Hz | 1.053 Hz (taut-string) | **+0.28%** |
| Tension target | 1.714 MN | 1.714 MN at t=0 | +0.008% |
| Dynamic tension range | - | 2.11 kN | - |
| Position response | - | mid-node disp. 0.0114 m, sag range 0.0114 m | - |

**整合性**: C1E は短い locked-coil cable (98.3 m) なので，taut-string から
推定した f1 は published f1 とほぼ一致する．入力諸元の整合性は良好．

**対比での価値**: この事例は Fujino universal curve が **有効な領域**の anchor として使える．
ただし張力は `T=(2Lf1)^2 m` から作った入力ターゲットなので，solver が
独立に文献張力を再同定したわけではない．
