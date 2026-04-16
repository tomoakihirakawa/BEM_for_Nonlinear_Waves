# NRC Canada Large-Scale RWIV Test

Primary role: controlled parametric rain-wind-induced vibration benchmark.

## Source

| Ref | Citation | DOI/URL | Access status |
|---|---|---|---|
| S1 | A. D'Auteuil, S. McTavish, A. Raeesi, G. Larose, "An investigation of rain-wind induced vibrations on stay cables in a novel range of operating conditions", Journal of Wind Engineering and Industrial Aerodynamics, 242, 105581, 2023 | https://doi.org/10.1016/j.jweia.2023.105581 | Open access; raw data not shared in article |

## Test Rig and Cable Model

| Quantity | Value | Unit | Source | Status |
|---|---:|---|---|---|
| Facility | NRC Canada large-scale RWIV dynamic test rig | - | S1 | published |
| Cable diameter | 0.218 | m | S1 abstract | published |
| Moving mass | 217 | kg | S1 abstract | published |
| Heave frequency | 1.23 | Hz | S1 abstract | published |
| Structural damping ratio | 0.05-0.11 | percent | S1 abstract | published |
| Tested cable inclinations | 30, 55 | deg | S1 abstract | published |
| Positive yaw range | up to 45 | deg | S1 abstract | published |
| Negative yaw range | down to -40 | deg | S1 abstract | published |
| RWIV wind speed range | 5-12 | m/s | S1 abstract | published |
| Maximum amplitude | 150 | mm | S1 abstract | published |

## Representative Input Values

The physical test used a rig rather than a full stay-cable span. The local
input uses a 12 m representative span so that the published moving mass maps
to a line density of 18.083 kg/m.

| Quantity | Value | Unit | Source | Status |
|---|---:|---|---|---|
| Representative span | 12 | m | local model | assumed |
| Cable inclination | 30 | deg | S1 | published test condition |
| Diameter | 0.218 | m | S1 | published |
| Line density | 18.083 | kg/m | local model | inferred from 217 kg / 12 m |
| Heave frequency | 1.23 | Hz | S1 | published |
| Mean tension | 1.576e4 | N | local model | inferred from `T=(2 L f)^2 m` |
| EA | 7.465e9 | N | local model | inferred from circular area and E=200 GPa |

## Files

| File | Purpose |
|---|---|
| `settings.json` | Standalone GUI/C++ input using the 30 deg representative test case |
| `cables/NRC_RWIV_30deg.json` | Per-cable detail file for the representative test span |

## Modeling Notes

- The rig boundary conditions are not a real bridge-cable boundary condition.
- `settings.json` distributes the mean tension into top and bottom endpoint
  tensions so the gravity-loaded static initial condition is internally
  consistent.
- The input is intended for controlled parametric checks after adding RWIV
  forcing; it is not a direct reproduction of the apparatus.
- Use the article figures for onset and amplitude digitization. Raw data are
  not available from the article.
