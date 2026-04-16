# Calculation Snapshots

Generated on 2026-04-16 from `/tmp/cable_verify_*` solver outputs.

The figures are SVG files made directly from solver `paraview/*.vtp` node
positions. They can be opened in a browser.

## 揺れが分かりやすい図

| Figure | What to look at |
|---|---|
| [`sutong_2008_pos_neg_neutral_x10.svg`](sutong_2008_pos_neg_neutral_x10.svg) | Sutong 3 cables. Dotted black = neutral/static shape, red = positive extreme, blue = negative extreme. Displacement from neutral is amplified 10x. |
| [`sutong_2008_C_547m_pos_neg_neutral_x10.svg`](sutong_2008_C_547m_pos_neg_neutral_x10.svg) | Clearest single cable. C_547m positive/negative extremes around the neutral shape, displacement amplified 10x. |
| [`sutong_2008_C_547m_pos_neg_neutral.svg`](sutong_2008_C_547m_pos_neg_neutral.svg) | Same C_547m plot at actual scale. Use this to avoid over-reading the amplified figure. |
| [`visible_response_ranking.svg`](visible_response_ranking.svg) | Whole-case ranking by direct mid-node sag range. |

## Published 12-Stay Validation

| Figure | What to look at |
|---|---|
| [`fred_hartman_table15_published_properties.svg`](fred_hartman_table15_published_properties.svg) | FHWA Table 15 source properties: length, tension, mass, frequency. |
| [`fred_hartman_table15_12stays_small_multiples.svg`](fred_hartman_table15_12stays_small_multiples.svg) | 12 stays shown separately so the lines do not collapse into one overplotted chord. |
| [`fred_hartman_table15_tension_error.svg`](fred_hartman_table15_tension_error.svg) | Error between equilibrium solver `T_top` and published Table 15 tension. |
| [`fred_hartman_table15_frequency_consistency.svg`](fred_hartman_table15_frequency_consistency.svg) | Error between published frequency and `f=(1/2L)sqrt(T/m)` from the same table. |
| [`fred_hartman_table15_summary.svg`](fred_hartman_table15_summary.svg) | Compact length/frequency trend summary for the 12 stays. |

## 誤差図

| Figure | Reference value |
|---|---|
| [`validation_tension_error.svg`](validation_tension_error.svg) | Error between solver `top_tension[0]` and each `cables/*.json` `tension_top` reference. |
| [`validation_frequency_error.svg`](validation_frequency_error.svg) | Error between input-derived taut-string `f1` and published/reference `f1` where available. |

## 全体まとめ

| Rank | Case | Cable | Direct sag range [m] | Tension target error [%] | Note |
|---:|---|---|---:|---:|---|
| 1 | sutong_2008 | C_547m | 0.6117 | +1.478 | Most visible current public demo |
| 2 | sutong_2008 | C_454m | 0.3358 | +1.493 | Visible, but smaller than C_547m |
| 3 | sutong_2008 | C_334m | 0.1224 | +1.364 | Tower period was tuned from this cable's nominal 2:1 condition |
| 4 | stonecutters_2009 | C_long | 0.0457 | -0.857 | Placeholder case; not paper-validated |
| 5 | east_china_sea_viv_2019 | CAC20 | 0.0407 | +0.005 | Geometry/frequency anchor; VIV lock-in not modeled |
| 6 | stonecutters_2009 | C_mid | 0.0242 | +0.296 | Placeholder case |
| 7 | dongting_lake_2007 | A12_representative | 0.0239 | +0.000 | Representative seed with assumptions |
| 8 | fred_hartman_1995 | C_long | 0.0233 | +0.044 | Placeholder case |
| 9 | stavanger_city_2021 | C1E | 0.0114 | +0.008 | Good geometry/frequency consistency |
| 10 | stonecutters_2009 | C_short | 0.0037 | -0.312 | Placeholder case |
| 11 | nrc_rwiv_2023 | NRC_RWIV_30deg | 0.0020 | +0.001 | Rig mapped to simplified cable span |
| 12 | fhwa_dry_inclined_2007 | FHWA_Dry_1B | 0.0017 | +0.002 | Aerodynamic instability not modeled |

The displacement metric is based on actual output node positions, not on
tension-to-amplitude back-calculation. The amplified plots multiply only the
dynamic displacement from the neutral/static shape; the actual-scale figure is
included to keep the magnitude honest.
