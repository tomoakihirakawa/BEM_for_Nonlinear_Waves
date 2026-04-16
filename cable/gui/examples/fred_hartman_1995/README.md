# Fred Hartman Bridge (Texas, USA, 1995)

米国テキサス州 Baytown の **2 経間ツインデッキ斜張橋**．
1995 年開通直後から **rain-wind induced vibration (RWIV)** が深刻で，
世界的に RWIV 対策研究のベンチマークとなった歴史的事例．

## 諸元 (主要公開値)

| 項目 | 値 | 出典 |
|---|---|---|
| 形式 | 2 経間ツインデッキ斜張橋，H 型主塔 2 基 | FHWA レポート |
| 中央径間 | 381 m | 同 |
| 主塔高 | 134 m (デッキから) | 同 |
| ケーブル本数 | **192 本** (48 × 2 プレーン × 2 デッキ) | 同 |
| ケーブル長範囲 | 40 – 180 m | 同 |
| ケーブル直径範囲 | 160 – 220 mm | 同 |
| 開通 | 1995 | 同 |
| 現地調査・測定期間 | 1997–2007 | Main, Jones, FHWA |

## 観測された振動

| 現象 | 振幅 (peak-to-peak) | 条件 | 出典 |
|---|---|---|---|
| RWIV | **最大 1.8 m** | 降雨，風速 10–15 m/s，風向 30–60° | Main & Jones 2001 |
| Vortex excitation | 0.3–0.5 m | 強風，乾燥 | Phelan 2006 |
| 長期 SHM | RMS 加速度 g 単位 | 長期，全気象条件 | Zuo & Jones 2010 |

**モード**: 1 次が主だが，一部ケーブルで 2–3 次モードも励振．
**頻度**: 開通 2 年で数十回の大振幅事象，damper 設置前は構造部材疲労懸念．

## ダンパー (retrofit)

1997 年以降，Infrastructure Engineering 製 **粘性ダンパー** を順次追加:

- 位置: anchor から `x_d / L ≈ 0.03` (最も短いケーブル以外)
- 設計減衰: 目標 ζ ≥ 0.5% (Fujino universal curve 基準)
- 実測改善: ζ = 0.1–0.2% → **0.8–2%** (cable ごとに異なる)

## データ取得予定 — 実施状況

- [x] FHWA-HRT-05-083 Appendix F Table 15 の 12 stays を JSON 化
  → `published_table15/`, `settings_table15_12stays.json`
- [ ] 全 192 本ケーブル諸元の表抽出可否を追加確認
- [ ] RWIV 時系列データ (Main & Jones 2001 Fig. 7–9) の数値化
- [ ] 気象データ (風速・雨量) 長期時系列
- [ ] ダンパー仕様 (type, c_d, x_d/L) テーブル化
- [ ] 数値再現: Fujino universal curve との比較

## 主要文献 (URL / DOI)

1. **Kumarasena, S., Jones, N.P., Irwin, P., Taylor, P. (2007)**
   *"Wind-Induced Vibration of Stay Cables"* **FHWA-HRT-05-083**.
   U.S. Department of Transportation, Federal Highway Administration. (約 300 頁)
   - Report page: https://www.fhwa.dot.gov/publications/research/infrastructure/structures/05083/
   - PDF archive (ROSA P): https://rosap.ntl.bts.gov/view/dot/50346
   - ★最重要・公開無料

2. **Main, J.A. & Jones, N.P. (2001)**
   "Full-scale measurements of stay cable vibration"
   *J. Sound Vib.* **245**(3), 435–456.
   - DOI: https://doi.org/10.1006/jsvi.2001.3555

3. **Zuo, D. & Jones, N.P. (2010)**
   "Interpretation of field observations of wind- and rain-wind-induced
   stay cable vibrations" *J. Wind Eng. Ind. Aerodyn.* **98**(2), 73–87.
   - DOI: https://doi.org/10.1016/j.jweia.2009.09.003

4. **Phelan, R.S., Sarkar, P.P., Mehta, K.C. (2006)**
   "Full-scale measurements to investigate rain-wind induced cable-stay
   vibrations" *J. Struct. Eng.* **132**(12), 1947–1956.
   - DOI: https://doi.org/10.1061/(ASCE)0733-9445(2006)132:12(1947)

## 研究テーマとの関係

- **Theme A (support-excited vs direct damper)**: Fred Hartman は
  **直接風荷重 (RWIV)** による振動で，塔・桁の連成は二次的．
  本研究の **塔経由支点励振 vs 直接風** の対比の「**直接風側の典型**」として使える．
- **Fujino universal curve の検証**: 192 本ケーブルの damper 追加前後の
  ζ 変化を universal curve で予測 → 実測と比較 → 破綻領域の同定．
- **由利橋の代替 anchor**: 社外秘制約なしでダンパー関連論文の背景に使える．

## Published 12-Stay Table 15 Model

FHWA-HRT-05-083 Appendix F, Table 15 gives 12 south-tower central-span stays
(`13S` to `24S`) with mass, tension, length, angular frequency, and frequency.
These values are now included in [`published_table15/`](published_table15/).

| File | Purpose |
|---|---|
| `published_table15/*.json` | One JSON per published Table 15 stay |
| `settings_table15_12stays.json` | Dynamic 12-stay settings using the existing tower placeholder |
| `settings_table15_12stays_equilibrium.json` | Equilibrium-only 12-stay target check |

Important limitation: Table 15 does not publish endpoint coordinates, cable
outer diameter, or EA. The JSON files therefore preserve the published
`mass/T/L/f` values, while endpoint geometry is an assumed 30 degree planar
chord and EA/diameter are inferred from steel-equivalent area.

Internet check (2026-04-16): the official FHWA HTML/PDF pages show Figure 117
for the 13S-24S three-dimensional arrangement and Table 15 for the published
properties, but the table itself contains no endpoint coordinates. A later
Figure 126 gives an equivalent 2D model with AS1-AS12 labels and restrainer
locations; that is a different A-line side-span data set from the 13S-24S
Table 15 set, so it should not be mixed into these 12 Table 15 JSON files.

### Verification (2026-04-16)

The table frequencies are internally consistent with
`f = (1 / (2 L)) sqrt(T / m)`. The maximum mismatch is 0.126%, consistent with
rounding in the published table.

Equilibrium and dynamic solver runs completed for all 12 stays:

| Check | Result |
|---|---:|
| Equilibrium convergence | 12/12 |
| Dynamic completion | 12/12 |
| Max equilibrium `T_top` error vs published T target | 0.819% |
| NaN/Inf/negative tension | 0 |
| Visible dynamic motion with current generic tower setting | negligible |

This makes Fred Hartman Table 15 a good 12-cable published-properties
validation case. It is not yet a calibrated vibration-amplitude reproduction.

### Figures

The following figures are generated from the computation results and published
Table 15 values.

![Fred Hartman Table 15 published properties](../../../docs/calculation_snapshots/fred_hartman_table15_published_properties.svg)

![Fred Hartman Table 15 small multiples](../../../docs/calculation_snapshots/fred_hartman_table15_12stays_small_multiples.svg)

![Fred Hartman Table 15 tension error](../../../docs/calculation_snapshots/fred_hartman_table15_tension_error.svg)

![Fred Hartman Table 15 frequency consistency](../../../docs/calculation_snapshots/fred_hartman_table15_frequency_consistency.svg)
