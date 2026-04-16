# Stonecutters Bridge (香港昂船洲大橋, Hong Kong, 2009)

香港の **1018 m 主径間斜張橋**．224 本ケーブル．
香港理工大 (PolyU) + 香港科技大 (HKUST) の **WASHMS** (Wind and Structural
Health Monitoring System) が世界でも有数の詳細な長期モニタリングを運用中．

## 諸元 (主要公開値)

| 項目 | 値 | 出典 |
|---|---|---|
| 形式 | 2 経間連続斜張橋，単柱型主塔 2 基 | Wong 2007 |
| 主径間 | 1018 m | 同 |
| 側径間 | 289 m (各) | 同 |
| 主塔高 | 298 m (海面から) | 同 |
| ケーブル本数 | **224 本** (4 × 56) | Ni et al. 2010 |
| ケーブル長範囲 | 82 – 540 m | 同 |
| ケーブル直径 | 140–160 mm | 同 |
| 開通 | 2009 | 同 |
| WASHMS 運用開始 | 2009 以降継続 | Wong 2007 |

## 観測された振動 / 長期モニタリング

| 項目 | 内容 | 出典 |
|---|---|---|
| 加速度時系列 | 全ケーブル，1–200 Hz サンプリング，年単位 | Ni 2011+ |
| 風速・方向 | デッキ・主塔複数点で 1 Hz 以上 | Wong 2007 |
| 温度 | ケーブル・主塔・deck 複数点 | Ko & Ni 2005 |
| **張力と温度の相関** | 年間変動 ±5% 程度 | Wang et al. 2014 |
| 減衰比の operational 同定 | モードごと，長期変動 | Sun & Ni 2017 |
| 大振幅事象 | 台風通過時の記録あり | Wong 2007 |

**特徴**: 振幅自体は Fred Hartman, Sutong ほど目立たないが，
**長期・連続・多パラメータ同時計測**という唯一無二のデータ資産．
**温度と damper performance の関係**を実測で追える．

## ダンパー

設計段階から **external viscous damper** を全 224 本に搭載:

- 位置: anchor から `x_d / L ≈ 0.03`
- 種類: 粘性ダンパー (各ケーブル最適化済み)
- 実測 ζ: 設計 0.5–1% を達成 (Sun & Ni 2017)

## データ取得予定 — 実施状況

- [ ] Wong (2007) WASHMS 公式レポート取得
- [ ] Ni, Ye, Ko (2010) *J. Eng. Mech.* 138 取得
- [ ] Sun & Ni (2017) operational modal ID 論文取得
- [ ] Wang et al. (2014) 温度効果論文取得
- [ ] 代表ケーブル (最長 ~540 m および中間 ~200 m) 諸元 JSON 化
- [ ] 温度-減衰相関データ抽出

## 主要文献 (URL / DOI)

1. **Wong, K.Y. (2007)**
   "Design of a structural health monitoring system for long-span bridges"
   *Structure and Infrastructure Engineering* **3**(2), 169–185.
   - DOI: https://doi.org/10.1080/15732470600591117
   - (WASHMS 公式設計書，Taylor & Francis)

2. **Ni, Y.Q., Ye, X.W., Ko, J.M. (2010)**
   "Monitoring-based fatigue reliability assessment of steel bridges:
   Analytical model and application"
   *J. Struct. Eng.* **136**(12), 1563–1573.
   - DOI: https://doi.org/10.1061/(ASCE)ST.1943-541X.0000250
   - Note: 論文タイトル/巻号は ASCE サイトで要確認．上記は著者グループの
     代表的な長期 SHM データ論文．

3. **Ko, J.M. & Ni, Y.Q. (2005)**
   "Technology developments in structural health monitoring of large-scale
   bridges" *Engineering Structures* **27**(12), 1715–1725.
   - DOI: https://doi.org/10.1016/j.engstruct.2005.02.020

4. **Sun, L. & Ni, Y.Q. (以降の operational modal ID 論文群)**
   - Google Scholar 検索: https://scholar.google.com/scholar?q=Sun+Ni+Stonecutters+operational+modal+identification+stay+cable

5. **Wang, Z.H., Xu, Y.L., Xia, Y. (2014)**
   "Temperature effects on vibrations of stay cables" *J. Sound Vib.* **333**(7), 1863–1880.
   - DOI: https://doi.org/10.1016/j.jsv.2013.11.029

## 研究テーマとの関係

- **Theme B (温度依存・低温性能劣化)**: 香港は年間温度レンジ 10–35 ℃ と
  狭いが，**温度-張力-減衰の三者関係の公開データ**が最も整っている．
  由利橋の温度依存仮説を **公開 anchor で補強**できる唯一の事例．
- **長期 validation データ**: 本ソルバでの年単位シミュレーション結果を
  実測長期時系列と比較 → モデル妥当性の強力な論拠．
- **多目的 SHM 同定**: 減衰比の変動を operational modal analysis で同定
  している研究がある → 減衰設計式の逆問題的検証に使える．
