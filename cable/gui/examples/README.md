# Examples

このディレクトリは `cable_solver` と pycable GUI で読める入力例を置く場所です．
この README は公開・共有可能な例の整理用で，`yuri_bridge/` は対象外です．

## 使い分け

| 種別 | 対象 | 用途 |
|---|---|---|
| 最小入力例 | `synthetic/` | ソルバ入力形式，GUI 表示，静的/動的モードの動作確認 |
| 公開文献ベンチマーク | `dongting_lake_2007/`, `stavanger_city_2021/`, `east_china_sea_viv_2019/`, `nrc_rwiv_2023/`, `fhwa_dry_inclined_2007/` | 論文・公的レポート値に基づく RWIV/VIV/乾燥傾斜ケーブル検証 seed |
| 長大橋ケースノート | `fred_hartman_1995/`, `stonecutters_2009/`, `sutong_2008/` | ダンパー，長期 SHM，多モード振動の公開情報整理と将来の詳細 JSON 化 |
| 作業テンプレート | `literature_benchmarks/` | 新規公開文献ケースの抽出項目，仮定，ファイル規約 |

`bridge_C01.json` は単一ケーブルのローカル回帰用 seed です．公開文献ベンチマークとしての出典整理が未完了なので，下の公開ケース表には入れていません．

## 最小入力例

`synthetic/` に置いている 4 つは公開文献の再現ケースではなく，ソルバ入力形式と GUI 動作を確認するための synthetic / smoke-test 用 seed です．
`synthetic/catenary_500m.json` と `synthetic/bridge_cable_small.json` は GUI テストからも参照されています．

| File | 入力形式 | 内容 |
|---|---|---|
| `synthetic/catenary_500m.json` | legacy single-line | 500 m 水平スパン，58 m 水深，重いチェーンの静的カテナリ．`cable_solver.cpp` の既定値に近い回帰確認用． |
| `synthetic/bridge_cable_small.json` | legacy single-line | 50 m 級の小さな橋梁ケーブル．橋梁スケールの静的つり合いを軽く確認する synthetic case． |
| `synthetic/dynamic_heave.json` | per-cable dynamic | 500 m カテナリの片端を鉛直加振する動的入力例．時刻歴，スナップショット，GUI 再生の確認用． |
| `synthetic/mooring_3leg_catenary.json` | multi-line BEM-compatible | `mooring_L01` から `mooring_L03` までの 3 本係留を 13 要素配列で与える BEM 互換形式の確認用． |
| `bridge_C01.json` | per-cable equilibrium | 橋梁ケーブル単体の回帰 seed．公開出典・代表性は未整理． |

## 公開文献ベンチマーク

各ディレクトリには，原則として `README.md`, `settings.json`, `cables/*.json` があります．
`settings.json` は GUI / C++ にそのまま渡せる実行単位，`cables/*.json` は個別ケーブルの値を確認するための詳細ファイルです．

| Directory | 主対象 | 代表入力 | 主な公開値 | ローカル推定・仮定 |
|---|---|---|---|---|
| `dongting_lake_2007/` | Dongting Lake Bridge の実橋 RWIV | `A12_representative` | 45 日計測，風速 6-14 m/s，相対 yaw 10-50 deg，雨量 < 8 mm/h，最大 10 g，0.7 m peak-to-peak，主に 3 次モード | 代表ケーブルの一部寸法，質量，固有振動数，張力を暫定設定 |
| `stavanger_city_2021/` | Stavanger City Bridge の乾湿比較 RWIV | `C1E` | C1E/C1W, C2E/C2W の長さ・径・質量・1 次振動数，湿潤時の 3 次モード大振幅 | 張力を `T=(2 L f_1)^2 m` から推定，EA を E=200 GPa 相当で推定 |
| `east_china_sea_viv_2019/` | 長大橋ケーブルの高次・多モード VIV | `CAC20` | CAC20 長さ 323 m，径 0.12 m，傾斜 26.3 deg，1 次 0.378 Hz，18-40 次モード | 質量 80 kg/m を代表値として仮定，張力・EA を推定 |
| `nrc_rwiv_2023/` | NRC 大型 RWIV 試験 | `NRC_RWIV_30deg` | 径 0.218 m，移動質量 217 kg，heave 1.23 Hz，傾斜 30/55 deg，yaw -40 から 45 deg，RWIV 5-12 m/s | 実験装置を 12 m 代表スパンに写像し，線密度を 217 kg / 12 m として設定 |
| `fhwa_dry_inclined_2007/` | FHWA/NRC 乾燥傾斜ケーブル風洞試験 | `FHWA_Dry_1B` | 長さ 6.7 m，直接曝露長 5.9 m 以上，径 160 mm，質量 60.8 kg/m，複数 damping level と yaw/inclination matrix | 1 Hz 代表振動数から張力を推定，EA を円断面・E=200 GPa で推定 |

## 長大橋ケースノート

次のディレクトリは公開文献・公的レポートに基づく調査メモと暫定 JSON を含みます．詳細値の表抽出や論文照合は継続中です．

| Directory | 主対象 | 現状 |
|---|---|---|
| `fred_hartman_1995/` | Fred Hartman Bridge の RWIV と粘性ダンパー retrofit | FHWA-HRT-05-083 Appendix F Table 15 の published 12 stays を `published_table15/` に追加済み．既存 short/mid/long は古い placeholder． |
| `stonecutters_2009/` | Stonecutters Bridge の WASHMS，温度-張力-減衰，長期 SHM | `tower.json` と短・中・長ケーブル seed あり．WASHMS 論文群の照合と温度・減衰データ抽出が未完了． |
| `sutong_2008/` | Sutong Bridge の長尺ケーブル，多モード VIV，viscous-shear damper | `tower.json` と 334/454/547 m ケーブル seed あり．多モードダンパー仕様と代表ケーブル表を抽出予定． |

## 公式・公開ソース

| Case | Primary source | Link | Access | この repo での扱い |
|---|---|---|---|---|
| FHWA/NRC dry inclined, Fred Hartman reference | Kumarasena, Jones, Irwin, Taylor, *Wind-Induced Vibration of Stay Cables*, FHWA-HRT-05-083, 2007 | https://www.fhwa.dot.gov/publications/research/infrastructure/bridge/05083/ | Public FHWA report | 乾燥傾斜試験，Scruton 数，RWIV 設計基準，Fred Hartman ダンパー情報の基準資料 |
| Dongting Lake Bridge | Ni, Wang, Chen, Ko, *Field observations of rain-wind-induced cable vibration in cable-stayed Dongting Lake Bridge*, JWEIA 95(5), 2007 | https://doi.org/10.1016/j.jweia.2006.07.001 | Abstract/metadata public; full text may need access | 実橋 RWIV の代表ケース．入力は現時点で calibration 済みではなく seed |
| Stavanger City Bridge | Daniotti, Jakobsen, Snæbjörnsson, Cheynet, Wang, *Observations of bridge stay cable vibrations in dry and wet conditions*, JSV 503, 2021 | https://doi.org/10.1016/j.jsv.2021.116106 | Open access | 湿潤/乾燥の比較と C1E 代表入力 |
| East China Sea long-span bridge | Chen, Gao, Laima, Li, *A Field Investigation on Vortex-Induced Vibrations of Stay Cables in a Cable-Stayed Bridge*, Applied Sciences 9(21), 4556, 2019 | https://doi.org/10.3390/app9214556 | Open access | 高次・多モード VIV の代表入力 |
| NRC large-scale RWIV | D'Auteuil, McTavish, Raeesi, Larose, *An investigation of rain-wind induced vibrations on stay cables in a novel range of operating conditions*, JWEIA 242, 105581, 2023 | https://doi.org/10.1016/j.jweia.2023.105581 | Open access | 制御された RWIV パラメータ試験の代表入力 |
| Stonecutters Bridge WASHMS | Wong, *Design of a structural health monitoring system for long-span bridges*, Structure and Infrastructure Engineering 3(2), 2007 | https://doi.org/10.1080/15732470600591117 | Publisher access | WASHMS と長期 SHM の入口資料．値の詳細抽出は未完了 |
| Sutong Bridge multimode damper | Chen et al., *Multimode cable vibration control using a viscous-shear damper: Case studies on the Sutong Bridge*, Structural Control and Health Monitoring 27(6), 2020 | https://doi.org/10.1002/stc.2536 | Publisher access | 多モード damping 設計の入口資料．ダンパー仕様抽出は未完了 |

## 出典ステータスの読み方

| Status | 意味 |
|---|---|
| `published` | 論文，FHWA レポート，公式ページに表・本文値として出ている値 |
| `inferred` | 公開値から `T=(2 L f)^2 m` などで計算した値 |
| `assumed` | ソルバ入力を閉じるために置いた暫定値．後で表抽出や同定で置換すべき値 |
| `local` | この repo 内の代表名，実行条件，出力先など，文献値ではない管理情報 |

公開ベンチマークとして使うときは，各ケースの README の `Status` 列を確認し，`assumed` や `local model` の値を実測再現の根拠として扱わないでください．

## 風・振動の設計式メモ

FHWA-HRT-05-083 は，斜張橋ケーブルの風起因振動を整理する基準資料として使えます．
この repo の JSON は構造ケーブル入力なので，以下は空力モデルや後処理を追加するときの参照式です．

| Quantity | Formula | 用途 |
|---|---|---|
| Reynolds number | `Re = rho V D / mu` | 風速，径，空気粘性から流れ領域を整理 |
| Strouhal number | `S = N_s D / V` | 渦励振周波数と風速の対応 |
| Scruton number | `Sc = m zeta / (rho D^2)` | 質量・減衰・径から風起因振動の抑制余裕を評価 |
| Vortex lock-in wind speed | `V = N_r D / S` | 固有振動数 `N_r` と Strouhal 数から渦励振風速を見積もる |
| RWIV mitigation criterion | `Sc > 10` | FHWA が示す雨風振動抑制の目安 |
| Wake/dry inclined criterion | `U_crit / (f D) = C sqrt(Sc)` | wake galloping / dry inclined galloping の概略安定判定．`C` は配置依存 |

これらの式は現在の `settings.json` に直接入っていません．現状の入力は，端点，長さ，線密度，EA，重力，減衰，風速ベクトルなどをソルバへ渡す構造モデルです．

## 入力形式

`cable_solver.cpp` は入力 JSON を自動判別します．

| Format | 判別キー | 代表ファイル | 説明 |
|---|---|---|---|
| settings mode | `input_files` | `fred_hartman_1995/settings.json` | 複数の cable/body JSON を読み込む実行単位．BEM 連成風の構成に近い． |
| per-cable | `end_a_position` | `synthetic/dynamic_heave.json`, `*/cables/*.json` | 1 本のケーブルを自己完結に定義．`end_a_body`, `end_b_body`, `initial_condition`, `tension_top` などを持てる． |
| multi-line BEM-compatible | `mooring_*` / `cable_*` | `synthetic/mooring_3leg_catenary.json` | BEM 入力に近い 13 要素配列形式． |
| legacy single-line | `point_a`, `point_b` | `synthetic/catenary_500m.json` | 旧形式の単線カテナリ入力． |

## 主なパラメータ

| Key | Unit | Meaning |
|---|---|---|
| `point_a`, `point_b` | m | legacy single-line の端点座標 |
| `end_a_position`, `end_b_position` | m | per-cable 形式の端点座標 |
| `cable_length` | m | 無伸長長さ．直線距離より短いと静的つり合いが成立しにくい |
| `n_segments`, `n_points` | - | 離散化数．`n_points = n_segments + 1` |
| `line_density` | kg/m | 線密度．水中係留では浮力補正後の実効値を使う場合がある |
| `EA` | N | 軸剛性 |
| `diameter` | m | 空力・流体力用の代表径 |
| `damping` | - | 現状ソルバの数値緩和/減衰係数．文献の構造減衰比とは同一視しない |
| `initial_condition` | string | `tension` なら `tension_top`, `tension_bottom` を初期条件調整に使う |
| `wind_type`, `wind_U_mean` | string, m/s | 一様風などの風条件 |
| `mode` | string | `equilibrium` または `dynamic` |
| `dt`, `t_end`, `output_interval` | s, s, step | 動的計算の時間刻み，終了時刻，出力間隔 |
| `max_equilibrium_steps`, `equilibrium_tol` | step, m/s | 静的つり合い計算の反復上限と収束判定 |
| `snapshot_interval` | step | GUI が読む `SNAPSHOT` 出力間隔 |

## 実行例

```bash
cd <cpp>/cable/build_solver
./cable_solver ../gui/examples/synthetic/catenary_500m.json /tmp/cable_out/
python -m json.tool /tmp/cable_out/result.json | head
```

pycable GUI からは `gui/run.sh` を起動し，`File -> Open...` で `settings.json` または単体 JSON を読み込みます．

```bash
cd <cpp>/cable/gui
./run.sh
```

## 追加時のルール

1. 公開文献ケースは，ケースディレクトリ直下に `README.md` を置き，出典，値，単位，`published`/`inferred`/`assumed` を明示する．
2. 文献値から計算した張力や EA は，式と前提を README に残す．
3. 生データや非公開案件は，この公開ベンチマーク索引に混ぜない．
4. 代表入力が calibration 済みでない場合は，`benchmark seed` と書き，再現済みケースのように扱わない．

---

# 独自再計算・検証結果 (2026-04-16)

以下は `build_solver/cable_solver` で 8 個の公開/公開候補ケースを再実行した結果です．
`yuri_bridge/` はこの公開ベンチマーク集計から除外しています．

```bash
for b in fred_hartman_1995 dongting_lake_2007 stavanger_city_2021 \
         east_china_sea_viv_2019 nrc_rwiv_2023 fhwa_dry_inclined_2007 \
         sutong_2008 stonecutters_2009; do
  out=/tmp/cable_verify_$b
  mkdir -p "$out"
  ./build_solver/cable_solver "gui/examples/$b/settings.json" "$out/" \
    > "/tmp/cable_verify_$b.log" 2>&1
done
```

## 結論

1. 全 8 ケース，14 ケーブルが完走し，`top_tension`, `bottom_tension`, `max_tension` に NaN/Inf/負値は出ませんでした．
2. `initial_condition: "tension"` のケースでは，出力の `top_tension(t=0)` は入力ターゲットにおおむね一致します．ただし，これは secant 法の初期条件合わせの確認であり，論文値との独立検証ではありません．
3. 実出力の節点位置から見た可視振動は Sutong が最大です．特に `C_547m` の中央 sag range が約 0.61 m で最も大きいです．
4. 以前の「張力変動から 2-3.5 m の中央振幅を逆算する」評価は過大です．現在の励振設定では，実位置出力から読む中央変位は Sutong でも 0.08-0.35 m 程度です．
5. 公開論文との比較で信頼できるのは，現時点では「入力諸元が文献値と整合しているか」と「動的計算が安定に完走するか」です．VIV/RWIV の実測振幅再現は，現ソルバに空力ロックイン/雨水リブレット/ギャロッピングモデルが未実装なので，まだ主張できません．

## 張力ターゲットとの一致

`Solver T0` は動的計算開始時の `top_tension[0]` です．`Target top` は各 `cables/*.json` に入っている `tension_top` です．

| Case | Cable | Target top [MN] | Solver T0 [MN] | Diff [%] | Tmax-Tmin [kN] |
|---|---|---:|---:|---:|---:|
| fred_hartman_1995 | C_long | 4.5000 | 4.5020 | +0.044 | 3.92 |
| fred_hartman_1995 | C_mid | 3.5000 | 3.4783 | -0.621 | 1.22 |
| fred_hartman_1995 | C_short | 2.5000 | 2.4867 | -0.531 | 0.29 |
| dongting_lake_2007 | A12_representative | 2.0358 | 2.0359 | +0.000 | 16.81 |
| stavanger_city_2021 | C1E | 1.7141 | 1.7142 | +0.008 | 2.11 |
| east_china_sea_viv_2019 | CAC20 | 4.8260 | 4.8263 | +0.005 | 9.94 |
| nrc_rwiv_2023 | NRC_RWIV_30deg | 0.0163 | 0.0163 | +0.001 | 2.05 |
| fhwa_dry_inclined_2007 | FHWA_Dry_1B | 0.0123 | 0.0123 | +0.002 | 1.61 |
| sutong_2008 | C_334m | 4.1825 | 4.2395 | +1.364 | 182.57 |
| sutong_2008 | C_454m | 5.0990 | 5.1751 | +1.493 | 200.50 |
| sutong_2008 | C_547m | 6.2405 | 6.3327 | +1.478 | 227.50 |
| stonecutters_2009 | C_long | 10.0000 | 9.9143 | -0.857 | 151.35 |
| stonecutters_2009 | C_mid | 5.5000 | 5.5163 | +0.296 | 98.22 |
| stonecutters_2009 | C_short | 3.0000 | 2.9906 | -0.312 | 27.07 |

Sutong の差が約 1.4-1.5% とやや大きいのは，入力コメントの `H` が水平張力成分として記録されている一方，solver の `top_tension` は端部セグメントの軸張力として出力されるためです．これは大きな破綻ではありませんが，`H` と `T_top` は同じ物理量として扱わない方が安全です．

## 実位置から見た可視振動ランキング

ここでは ParaView 出力または `SNAPSHOT` の節点位置から，中央節点の moving chord からの距離を計算し，その range を `sag range` として比較しました．張力変動から逆算した推定値ではありません．

| Rank | Case | Cable | Frames | Mid-node disp. [m] | Max-node disp. [m] | Sag range [m] |
|---:|---|---|---:|---:|---:|---:|
| 1 | sutong_2008 | C_547m | 21 | 0.3449 | 0.3449 | **0.6117** |
| 2 | sutong_2008 | C_454m | 21 | 0.1815 | 0.1818 | **0.3358** |
| 3 | sutong_2008 | C_334m | 21 | 0.0790 | 0.0790 | **0.1224** |
| 4 | stonecutters_2009 | C_long | 21 | 0.0457 | 0.0470 | 0.0457 |
| 5 | east_china_sea_viv_2019 | CAC20 | 41 | 0.0407 | 0.0407 | 0.0407 |
| 6 | stonecutters_2009 | C_mid | 21 | 0.0229 | 0.0245 | 0.0242 |
| 7 | dongting_lake_2007 | A12_representative | 41 | 0.0239 | 0.0240 | 0.0239 |
| 8 | fred_hartman_1995 | C_long | 21 | 0.0233 | 0.0233 | 0.0233 |
| 9 | stavanger_city_2021 | C1E | 41 | 0.0114 | 0.0114 | 0.0114 |
| 10 | fred_hartman_1995 | C_mid | 21 | 0.0086 | 0.0086 | 0.0088 |
| 11 | stonecutters_2009 | C_short | 21 | 0.0033 | 0.0033 | 0.0037 |
| 12 | fred_hartman_1995 | C_short | 21 | 0.0023 | 0.0023 | 0.0025 |
| 13 | nrc_rwiv_2023 | NRC_RWIV_30deg | 101 | 0.0020 | 0.0020 | 0.0020 |
| 14 | fhwa_dry_inclined_2007 | FHWA_Dry_1B | 101 | 0.0017 | 0.0044 | 0.0017 |

GUI で一番動きが見える候補は Sutong です．ただし，現在の `sutong_2008/tower.json` は「論文で観測された風応答を再現する設定」ではなく，塔の cantilever 変位を与えたデモ励振です．

## 周波数チェック

`top_tension` 時系列の FFT ピークを確認しました．これは厳密な固有振動数同定ではなく，現在の励振に対する応答周波数です．

| Case | Cable | Dominant tension peak [Hz] | Interpretation |
|---|---|---:|---|
| sutong_2008 | C_334m | 0.199 | tower period 5.08 s の強制成分が支配的 |
| sutong_2008 | C_454m | 0.199 | tower 成分が支配的．副ピーク 0.299 Hz は入力コメントの f1 と整合 |
| sutong_2008 | C_547m | 0.199 | tower 成分が支配的．副ピーク 0.249 Hz は入力コメントの f1=0.258 Hz に近い |
| stonecutters_2009 | C_long/C_mid/C_short | 0.995 | tower period 1.0 s の強制成分 |
| stavanger_city_2021 | C1E | 1.047 | 公開 f1=1.05 Hz と整合．ただし張力ターゲットはこの f1 から推定 |
| east_china_sea_viv_2019 | CAC20 | 0.399 | 公開 f1=0.378 Hz に近い．ただし質量 80 kg/m は仮定 |
| dongting_lake_2007 | A12_representative | 0.698 | local model の f1=0.65 Hz に近い．代表値は仮定を含む |

以前の「張力 FFT ピークを半分にして f1 とみなす」処理は，このデータには適用しません．Sutong/Stonecutters では塔の強制周期が支配的であり，自由振動同定ではありません．

## 文献整合性と限界

| Case | 文献値との関係 | 現時点の評価 |
|---|---|---|
| `sutong_2008/` | `cables/*.json` は Chen et al. 2020 Table 1 とされる L, m, D, H, A, inclination, λ² を入力化 | 形状・張力入力の整合チェックは可能．ただし H と solver axial tension の比較は物理量が完全には一致しない．動的振幅はデモ励振で，論文再現ではない |
| `stavanger_city_2021/` | Daniotti et al. 2021 Table 1 の L, D, m, f1 を利用 | taut-string から逆算した張力と f1 は整合．ただし張力は公開直接値ではなく入力推定値 |
| `east_china_sea_viv_2019/` | Chen et al. 2019 の CAC20 L=323 m, D=0.12 m, inclination=26.3 deg, f1=0.378 Hz を利用 | f1 と入力張力は整合．ただし mass per unit length は仮定．VIV ロックイン再現ではない |
| `dongting_lake_2007/` | Ni et al. 2007 の event range を代表化 | 代表ケーブルの寸法・質量・f1 は仮定が多く，定量 validation には未使用 |
| `nrc_rwiv_2023/` | NRC 大型 RWIV 試験値を 12 m span に写像 | 実験装置境界条件と異なるため，現状は入力 seed |
| `fhwa_dry_inclined_2007/` | FHWA/NRC 乾燥傾斜試験値を代表化 | 1 Hz 仮定で張力を置いた seed．Scruton 数・空力不安定の検証は未実装 |
| `fred_hartman_1995/` | FHWA-HRT-05-083 Appendix F Table 15 の 12 stays は mass, tension, length, omega, frequency が published | 12 本の published-properties validation が可能．ただし端点座標/EA/外径は未掲載なので geometry は仮定 |
| `stonecutters_2009/` | WASHMS/長期 SHM のケースノート | 現 JSON は placeholder．論文値との定量 validation には未使用 |

## 今回の発見

1. ソルバは，現在の公開ケース入力に対して安定に完走し，張力出力の数値破綻はありません．
2. ただし，`tension_top` を文献値または推定値として入力しているため，`T0` が文献値に近いことをもって「張力を独立再現した」とは言えません．
3. 可視振動の最大ケースは Sutong で変わりませんが，実位置ベースの振幅は張力逆算値よりかなり小さく，最大でも `C_547m` の sag range 約 0.61 m です．
4. 文献と整合していると言える強い点は，Stavanger と East China Sea の公開 f1/形状から作った入力が，taut-string 周波数と矛盾していないことです．
5. 研究用 validation として次に必要なのは，自由振動テストで f1 を直接同定し，文献 f1 と比較することです．現在の強制振動結果は tower/wind 入力に支配されます．

## 推奨追加検証

1. **自由振動 f1 同定**: 静的平衡後に中央節点へ小変位を与え，無風・無強制で解放し，変位 FFT で f1 を取る．
2. **文献振幅の再現は別タスク化**: VIV/RWIV には lock-in，rain rivulet，galloping などの空力不安定モデルが必要．
3. **Sutong の H/T 定義整理**: Chen 2020 Table 1 の `H` が水平張力なら，solver 側の axial `top_tension` と別列で扱う．
4. **placeholder ケースの格下げ表示**: Fred Hartman は Table 15 の published 12 stays と古い short/mid/long placeholder を分けて扱う．Stonecutters は詳細諸元抽出が終わるまで validation 済みケースとして扱わない．
