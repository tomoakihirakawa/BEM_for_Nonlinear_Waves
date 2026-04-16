# East China Sea Long-Span Bridge VIV Field Investigation

Primary role: high-mode and multimode vortex-induced vibration benchmark.

## Source

| Ref | Citation | DOI/URL | Access status |
|---|---|---|---|
| S1 | W.-L. Chen, D. Gao, S. Laima, H. Li, "A Field Investigation on Vortex-Induced Vibrations of Stay Cables in a Cable-Stayed Bridge", Applied Sciences, 9(21), 4556, 2019 | https://doi.org/10.3390/app9214556 | Open access |

## Bridge and Monitoring System

| Quantity | Value | Unit | Source | Status |
|---|---:|---|---|---|
| Bridge location | East China Sea | - | S1 | published |
| Main span | 620 | m | S1 | published |
| Deck width | 30.1 | m | S1 | published |
| Deck height | 3.0 | m | S1 | published |
| Stay cables | 168 | - | S1 | published |
| Monitored cables | 20 | - | S1 | published |
| Accelerometer sampling | 100 | Hz | S1 | published |
| Anemometer sampling | 32 | Hz | S1 | published |
| Analysis months | April, July, August 2010 | - | S1 | published |
| Validation months | September, October 2010 | - | S1 | published |

## Focus Cable CAC20

| Quantity | Value | Unit | Source | Status |
|---|---:|---|---|---|
| Cable length | 323 | m | S1 | published |
| Diameter | 0.12 | m | S1 | published |
| Aspect ratio | 2692 | - | S1 | published |
| Inclination | 26.3 | deg | S1 | published |
| Fundamental frequency | 0.378 | Hz | S1 | published |
| Mass per unit length | 80 | kg/m | local model | assumed typical stay value |
| Mean tension | 4.770e6 | N | local model | inferred from `T=(2 L f1)^2 m` |
| EA | 2.262e9 | N | local model | inferred from circular area and E=200 GPa |

## Observed VIV Response

| Quantity | Value | Unit | Source | Status |
|---|---:|---|---|---|
| Mean wind speed during large response | 4-6 | m/s | S1 | published |
| Broad large-response wind range | 3-12 | m/s | S1 | published |
| Reduced velocity during large response | about 5 | - | S1 | published |
| Main participative modes | 18-40 | - | S1 | published |
| Fundamental frequency used for modal index | 0.378 | Hz | S1 | published |
| Large response condition | low turbulence, wind nearly perpendicular to bridge axis | - | S1 | published |

## Files

| File | Purpose |
|---|---|
| `settings.json` | Standalone GUI/C++ input using CAC20 and 5.75 m/s uniform wind |
| `cables/CAC20.json` | Per-cable detail file using published CAC20 geometry and inferred mass/tension |

## Modeling Notes

- This is a dry-wind VIV case, not an RWIV case.
- `settings.json` distributes the mean tension into top and bottom endpoint
  tensions so the gravity-loaded static initial condition is internally
  consistent.
- The present solver input does not reproduce vortex lock-in by itself. It
  provides the structural cable model and the published wind scale for later
  addition of aerodynamic excitation models.
- The bridge name is intentionally kept as described in S1 until a primary
  source explicitly identifies it.

## Verification (2026-04-16)

Recomputed with:

```bash
./build_solver/cable_solver gui/examples/east_china_sea_viv_2019/settings.json /tmp/cable_verify_east_china_sea_viv_2019/
```

| Quantity | Published / input basis | Solver / derived value | 誤差 |
|---|---:|---:|---:|
| Cable length L | 323 m | 323 m | exact |
| Diameter D | 0.12 m | 0.12 m | exact |
| Inclination | 26.3° | 26.3° (end_a/end_b から復元) | exact |
| First natural frequency f₁ | 0.378 Hz | 0.380 Hz (taut-string, 仮定 m=80 kg/m) | **+0.59%** |
| Tension target | 4.826 MN | 4.826 MN at t=0 | +0.005% |
| Dynamic tension range | - | 9.94 kN | - |
| Position response | - | mid-node disp. 0.0407 m, sag range 0.0407 m | - |

**整合性**: CAC20 は L=323 m，傾斜 26° で張力が高く，taut-string から
推定した f1 は published f1 と近い．入力諸元の整合性は良好．

**注意**: mass per unit length m = 80 kg/m は S1 で未公開のため **inferred**．仮定値の妥当性は f₁ = 0.378 Hz との整合で逆算できるが，諸元表には `assumed` として記載．

**用途**: Theme C (multimode VIV) の open-access geometry/frequency anchor．
ただし現 solver 入力は quasi-steady wind であり，VIV lock-in や高次モード応答の
実測振幅を再現するものではない．
