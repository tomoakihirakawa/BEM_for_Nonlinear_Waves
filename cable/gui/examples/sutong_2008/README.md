# Sutong Bridge (苏通大桥, China, 2008)

中国江蘇省の **世界最長クラス斜張橋** (開通当時)．
主塔間 1088 m，ケーブル最大長 **577 m**．多モード VIV とマルチモード減衰
ダンパー設計の最新研究対象．

## 諸元 (主要公開値)

| 項目 | 値 | 出典 |
|---|---|---|
| 形式 | 2 経間連続斜張橋，倒 Y 型主塔 2 基 | Chen 2020 |
| 主径間 | 1088 m | 同 |
| 主塔高 | 300.4 m | 同 |
| ケーブル本数 | **272 本** (4 プレーン × 68) | 同 |
| ケーブル長範囲 | 84 – 577 m | Yan 2015 |
| ケーブル直径 | 140 – 170 mm | 同 |
| 開通 | 2008 | 同 |

## 観測された振動

| 現象 | 振幅 | 条件 | 出典 |
|---|---|---|---|
| 多モード VIV | 長尺ケーブルで低次 (n=1–5) 同時励振 | 風速 5–15 m/s | Zhou & Xu 2007 |
| 低周波ピーク振幅 | 1 次 + 3 次連成で最大 **RMS 0.2 m** | 長期 SHM | Yan 2015 |
| 風-温度相関 | 温度変動で張力変動，減衰特性変化 | 年間長期データ | Wang et al. 2012 |

**特徴**: 単一モード VIV ではなく **複数モードが同時に大振幅**．
Fujino の単モード universal curve では設計不可能 → multimode shear damper
を開発・実装した．

## ダンパー (設計段階で装備)

設計時から **MR fluid / viscous shear damper** を搭載:

- 位置: anchor から `x_d / L ≈ 0.03–0.04`
- 設計目標: 1–5 次モードに対して同時に $\zeta_n > 0.5\%$
- 実測: before damper ~0.15%, after damper **0.8–1.8%** (mode により差)

## データ取得予定 — 実施状況

- [ ] Chen, Sun, Huang (2020) *Struct. Control Health Monit.* 27(6), e2536 取得
- [ ] Yan 2015 SHM データ取得
- [ ] 272 本ケーブル諸元の代表値テーブル化 → `cables/` に最長 12 本のサンプル JSON 化
- [ ] multimode damper 仕様 (type, c_d(ω), x_d/L) 抽出
- [ ] 数値再現: multimode $\zeta_n$ の universal chart 拡張

## 主要文献 (URL / DOI)

1. **Chen, L., Sun, L., Xu, Y.L., Xu, F., Mai, X.Y. (2020)**
   "Multimode cable vibration control using a viscous-shear damper:
   Case studies on the Sutong Bridge"
   *Structural Control and Health Monitoring* **27**(6), e2536.
   - DOI: https://doi.org/10.1002/stc.2536
   - ★最重要★ (Wiley，大学アクセス要)

2. **Zhou, H.F., Xu, Y.L. (2007)**
   "Wind-induced vibrations of long stay cables"
   *J. Sound Vib.* **305**(1–2), 44–59.
   - DOI: https://doi.org/10.1016/j.jsv.2007.03.083

3. **Yan, B., Li, Q.Q. (2015)**
   "Long-term vibration monitoring of stay cables in Sutong Bridge"
   *Advances in Structural Engineering* **18**(11), 1821–1835.
   - DOI: https://doi.org/10.1260/1369-4332.18.11.1821

4. **Wang, Z.H., Xu, Y.L., Xia, Y. (2014)**
   "Temperature effects on vibrations of stay cables"
   *J. Sound Vib.* **333**(7), 1863–1880.
   - DOI: https://doi.org/10.1016/j.jsv.2013.11.029

## 研究テーマとの関係

- **Theme C (multimode hybrid damper chart)**: Sutong は **multimode 減衰
  設計が一次目標**の橋梁で，Fujino 単モード式の純粋な拡張として universal
  chart 拡張の対象に最適．
- **温度依存性データ**: Wang 2014 の温度-減衰相関データが取れれば
  **温度依存 Fujino 式** (Theme B) の検証にも使える．
- **長尺ケーブル**: 577 m の非常に長い cable の大変位挙動が実測されている
  → Krenk 漸近解の破綻領域確認に有利．

## Verification (2026-04-16)

Recomputed with:

```bash
./build_solver/cable_solver gui/examples/sutong_2008/settings.json /tmp/cable_verify_sutong_2008/
```

### Input consistency with Chen 2020 Table 1

The three cable JSON files encode the Table 1 values listed in their comments:
`L`, `m`, `D`, `H`, cross-sectional area, inclination, and `lambda^2`.
The solver uses `tension_top`/`tension_bottom` as initial-condition targets, so
the dynamic output at `t=0` is a target-consistency check, not an independent
identification of published tension.

| Cable | Input H/target top [kN] | Solver T_top(t=0) [kN] | Diff [%] | Note |
|---|---:|---:|---:|---|
| C_334m | 4,182.5 | 4,239.5 | +1.36 | solver output is axial end-segment tension |
| C_454m | 5,099.0 | 5,175.1 | +1.49 | same caveat |
| C_547m | 6,240.5 | 6,332.7 | +1.48 | same caveat |

The small positive offset is consistent with comparing a published horizontal
tension-like quantity `H` against the solver's axial segment tension. Do not
treat this table as a full paper validation.

### Dynamic demo response

`tower.json` imposes a cantilever motion with period 5.08 s. This is a visual
demo/support-excitation case, not a calibrated reproduction of the reported
wind-induced vibration in Chen 2020.

The table below uses actual output positions, not a back-calculation from
tension variation. `Sag range` is the range of the mid-node distance from the
instantaneous chord between the two endpoints.

| Cable | T0 [MN] | Tmax-Tmin [kN] | Mid-node disp. [m] | Max-node disp. [m] | Sag range [m] |
|---|---:|---:|---:|---:|---:|
| C_334m | 4.2395 | 182.57 | 0.0790 | 0.0790 | 0.1224 |
| C_454m | 5.1751 | 200.50 | 0.1815 | 0.1818 | 0.3358 |
| C_547m | 6.3327 | 227.50 | 0.3449 | 0.3449 | 0.6117 |

Sutong remains the most visible public demo case in this repository, but the
direct position-based motion is about 0.1-0.6 m in sag range, not the 2-3.5 m
that results from applying a taut-string tension/amplitude back-calculation to
this support-excited response.

### Frequency check

The dominant `top_tension` FFT peak is 0.199 Hz for all three cables, matching
the imposed tower period (5.08 s). Secondary peaks near published cable
frequencies appear for C_454m and C_547m, but this run is forced response, not
free-vibration identification. A proper validation of the published `f1` values
requires a separate free-vibration test.
