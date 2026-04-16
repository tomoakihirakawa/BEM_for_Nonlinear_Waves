# Input JSON Schema

`cable_solver` は 4 種類の入力形式を自動判別する．典型用途は
**per-cable 形式**（1 ケーブル = 1 JSON）または，複数ケーブルを束ねる
**settings.json 形式**である．

| セクション | 用途 |
|---|---|
| [Per-Cable 形式](#per-cable-形式推奨) | 単一ケーブルの完全自己完結 JSON |
| [settings.json](#settingsjson複数ケーブルの一括指定) | 複数 per-cable JSON + 共有スカラー |
| [動的モード](#動的モードprescribed-motion) | 端点強制振動による時間応答 |
| [流体と風](#流体と風fluid--wind) | 水／空気切替，一様・AR(1) 乱流風 |
| [単一ライン形式](#単一ライン形式レガシー) | レガシー フラットキー |
| [多重ライン形式](#多重ライン形式bem-互換) | BEM `mooring_<name>` 13 要素配列 |

## Per-Cable 形式（推奨）

1 ケーブル = 1 JSON ファイル．ケーブルが自分の端点位置，材料特性，初期条件を自己完結で持つ．

```json
{
  "name": "C01",
  "end_a_position": [-104.1, 0.0, 0.0],
  "end_b_position": [0.0, 0.0, 45.0],
  "end_a_body": "deck",
  "end_b_body": "tower",
  "initial_condition": "tension",
  "tension_top": 1019998,
  "tension_bottom": 986767,
  "cable_length": 113.41,
  "n_points": 31,
  "line_density": 63.2514,
  "EA": 1530925000.0,
  "damping": 0.5,
  "diameter": 0.10129,
  "output_directory": "./output",
  "gravity": 9.81,
  "mode": "equilibrium"
}
```

| Key | 必須 | Default | 意味 |
|---|---|---|---|
| `name` | いいえ | ファイル名 stem | ケーブル名 |
| `end_a_position` | はい | - | 端点 A の world 位置 [x, y, z] |
| `end_b_position` | はい | - | 端点 B の world 位置 [x, y, z] |
| `end_a_body` | いいえ | - | BEM 連携用 body ラベル |
| `end_b_body` | いいえ | - | BEM 連携用 body ラベル |
| `initial_condition` | いいえ | `"length"` | `"length"`: `cable_length` = 自然長．`"tension"`: 目標張力から自然長を反復で逆算 |
| `tension_top` | tension 時 | - | 上端の目標張力 [N] |
| `tension_bottom` | tension 時 | - | 下端の目標張力 [N] |
| `cable_length` | はい | - | `"length"` 時は自然長，`"tension"` 時はコード長（参照用） |
| `n_points` | はい | - | 節点数 (= `n_segments + 1`) |
| `line_density` | はい | - | kg/m |
| `EA` | はい | - | 軸剛性 [N] |
| `damping` | いいえ | 0.5 | N s/m |
| `diameter` | いいえ | 0.132 | [m] |
| `output_directory` | いいえ | - | 結果出力先（GUI が読み取り） |
| `gravity` | いいえ | 9.81 | m/s² |
| `mode` | いいえ | `"equilibrium"` | `"equilibrium"` / `"dynamic"` |
| `end_a_motion` | いいえ | `"fixed"` | 端点の運動（dynamic 用） |
| `end_b_motion` | いいえ | `"fixed"` | `"fixed"` / `"sinusoidal"` / `"cantilever"` |
| `end_b_motion_dof` | sinusoidal 時 | - | `"surge"` / `"sway"` / `"heave"` |
| `end_b_motion_amplitude` | sinusoidal 時 | - | 振幅 [m] |
| `end_b_motion_frequency` | sinusoidal 時 | - | 周波数 [Hz] |
| `dt` | dynamic 時 | - | 時間刻み [s] |
| `t_end` | dynamic 時 | - | 終了時刻 [s] |

## 初期条件: 自然長 vs 張力

- **`"initial_condition": "length"`**（デフォルト）: `cable_length` をそのまま自然長として使う．
- **`"initial_condition": "tension"`**: `tension_top` / `tension_bottom` を目標値として secant method で自然長を反復収束させる．両端張力の RMS 誤差を最小化する．実橋のように自然長が不明で張力のみ既知の場合に使う．

## settings.json（複数ケーブルの一括指定）

BEM の `settings.json` → `input_files` と同じパターン．

```json
{
  "input_files": ["cables/C01.json", "cables/C02.json", "cables/C12.json"],
  "output_directory": "./output",
  "gravity": 9.81,
  "mode": "equilibrium",
  "max_equilibrium_steps": 500000,
  "equilibrium_tol": 0.01,
  "snapshot_interval": 5000
}
```

- `input_files` の各パスは `settings.json` からの相対パスで解決される．
- settings の共通パラメータは各ケーブル JSON のデフォルトとして適用される．ケーブル側に同キーがあればそちらが優先される．
- 各ケーブルは独立に解かれ，結果は `output_dir/<name>_result.json` として個別生成される．

## 動的モード（prescribed motion）

端点に強制振動を与え，ケーブルの時間応答を解く．

```json
{
  "name": "D01",
  "end_a_position": [500.0, 0.0, -58.0],
  "end_b_position": [0.0, 0.0, 0.0],
  "end_b_motion": "sinusoidal",
  "end_b_motion_dof": "heave",
  "end_b_motion_amplitude": 2.0,
  "end_b_motion_frequency": 0.05,
  "cable_length": 522.0,
  "n_points": 41,
  "line_density": 348.5,
  "EA": 1400000000.0,
  "mode": "dynamic",
  "dt": 0.005,
  "t_end": 20.0,
  "output_interval": 0.1
}
```

動的結果 JSON には時系列の張力が含まれる．

```json
{
  "name": "D01",
  "mode": "dynamic",
  "n_output_steps": 201,
  "time": [0.0, 0.1, 0.2, "..."],
  "top_tension": ["T0", "T1", "..."],
  "bottom_tension": ["T0", "T1", "..."],
  "positions_final": [["x", "y", "z"], "..."],
  "computation_time_ms": 648.0
}
```

内部では CFL 安定条件
`dt_cfl = L0 / sqrt(EA / rho)` でサブステップを自動分割し，ユーザ指定の
`dt` が CFL を超えても安定に計算する．

## 流体と風（Fluid & Wind）

ケーブルが置かれる流体（水／空気）と，外部流速場を `settings.json` の
トップレベルで一括指定できる．system-wide 設定であり，per-cable 上書きはない．
全キー省略時は水中・無風で，レガシー入力と互換になる．

```json
{
  "input_files": ["cables/C01.json", "..."],
  "gravity": 9.81,
  "mode": "dynamic",
  "dt": 0.01,
  "t_end": 60.0,

  "fluid": "air",
  "fluid_density": 1.225,
  "drag_Cd": 1.2,

  "wind_type": "AR1",
  "wind_U_mean": [15.0, 0.0, 0.0],
  "wind_turbulence_intensity": 0.15,
  "wind_integral_time_scale": 5.0,
  "wind_seed": 12345
}
```

| Key | 値 | 意味 |
|---|---|---|
| `fluid` | `"water"` / `"air"` | プリセット選択．既定は `"water"` |
| `fluid_density` | number | 流体密度 [kg/m³]．プリセット上書き |
| `drag_Cd` | number | 抗力係数．プリセット上書き |
| `wind_type` | `"none"` / `"uniform"` / `"AR1"` | 風モデル．既定は `"none"` |
| `wind_U_mean` | `[Ux, Uy, Uz]` | 平均風速ベクトル [m/s] |
| `wind_turbulence_intensity` | number | `sigma_u / |U_mean|`．AR1 のみ |
| `wind_integral_time_scale` | number | `T_L` [s]．AR1 の OU 時定数 |
| `wind_seed` | integer | RNG 再現性．省略時は time-based |

**プリセット既定値**:

| プリセット | `fluid_density` [kg/m³] | `drag_Cd` [-] | 想定用途 |
|---|---|---|---|
| `"water"` | 1000.0 | 2.5 (Palm 2016) | 海洋係留ライン |
| `"air"` | 1.225 | 1.2 (亜臨界円柱) | 橋梁 stay，架空線 |

プリセット選択で `fluid_density` / `drag_Cd` の既定値が決まる．明示キーが与えられていれば，それがプリセットを上書きする．

**風モデル**（[WindField.hpp](../../lib/include/WindField.hpp)）:

| `wind_type` | 生成する風速 `U(X, t)` | 用途 |
|---|---|---|
| `"none"` | `0` | 無風（既定） |
| `"uniform"` | `U_mean`，時刻・位置ともに一定 | 静的たわみ，定常風下の動的応答 |
| `"AR1"` | `U_mean + x(t)`，`x` は Ornstein-Uhlenbeck 過程 | バフェッティング応答（準定常近似） |

AR(1) は全節点で同じ時系列を返す．つまり，現状は空間相関なしではなく，
完全相関の一様乱流である．1 本のケーブル応答を見る MVP としては妥当だが，
Davenport / von Kármán スペクトルや空間コヒーレンスは未実装である
（[風荷重について.md](../風荷重について.md) の Level 1 以降を参照）．

**適用の仕組み**:

抗力式は [Network.hpp](../../lib/include/Network.hpp) の
`getDragForce(Cd, rho, U_fluid)` を使う．

```text
F_d = 0.5 rho |v_rel|^2 Cd A_proj v_rel_hat
v_rel = U_fluid(X, t) - v_node
```

`cable_solver.cpp::resolveFluidConfig` が `settings.json` からプリセット解決と
風場クロージャ構築を行い，`applyFluidConfig` が全ケーブルの
`LumpedCable::FluidDensity`, `DragForceCoefficient`, `wind_field` メンバへ流し込む．

**平衡計算時の扱い**:

`solveEquilibrium()` は，風なし時は RAII ガードで `Cd=1000` / `rho=rho_water`
に差し替えて高速擬似緩和する．風あり時はこのガードを無効化し，物理的な
風下平衡に収束させる．動的モードでは，初回平衡を風を抑制して求め，
時間ループで風を再活性化する．

**典型的な使い方**:

- 海洋係留（水中，無風）: キー全省略，または `"fluid": "water"` のみ．
- 橋梁 stay の静的たわみ: `"fluid": "air"` + `"wind_type": "uniform"` + `"mode": "equilibrium"`．
- 橋梁 stay のバフェッティング応答: `"fluid": "air"` + `"wind_type": "AR1"` + `"mode": "dynamic"`．`t_end` は `T_L` より十分長く取る．

## 単一ライン形式（レガシー）

すべてのキーはオプションで，省略時はデフォルトが適用される．

| Key | Default | Unit | 備考 |
|---|---|---|---|
| `point_a` | `[500, 0, -58]` | m | 端点 A（アンカー想定） |
| `point_b` | `[0, 0, 0]` | m | 端点 B（フェアリード想定） |
| `cable_length` | `522` | m | 自然長の合計 |
| `n_segments` | `40` | - | 分割数．節点数は +1 |
| `line_density` | `348.5` | kg/m | 単位長あたり質量．旧引数名 `density` は同値 |
| `EA` | `1.4e9` | N | 軸剛性 |
| `damping` | `0.5` | N s/m | 隣接節点間相対速度に作用 |
| `diameter` | `0.132` | m | 抗力射影面積計算用 |
| `gravity` | `9.81` | m/s² | 正値で与える．向きは内部で -z |
| `mode` | `"equilibrium"` | - | 現状は `equilibrium` のみ受け付け |
| `max_equilibrium_steps` | `500000` | - | 上限ステップ数 |
| `equilibrium_tol` | `0.01` | m/s | 収束判定閾値（節点最大速度） |
| `snapshot_interval` | `10000` | - | SNAPSHOT 出力間隔 |

## 多重ライン形式（BEM 互換）

`mooring_<name>` または `cable_<name>` をキーに，13 要素のフラット配列を値として与える．
1 個でも存在すれば多重ラインモードに切り替わる．

```json
{
  "mooring_L01": [
    "L01",
    500.0, 0.0, -58.0,
    0.0, 0.0, 0.0,
    522.0,
    41,
    348.5,
    1.4e9,
    0.5,
    0.132
  ],
  "mooring_L02": ["L02", "..."],
  "gravity": 9.81,
  "mode": "equilibrium",
  "max_equilibrium_steps": 500000,
  "equilibrium_tol": 0.01,
  "snapshot_interval": 10000
}
```

配列要素の順序は
`[name, ax, ay, az, bx, by, bz, total_length, n_points, line_density, EA, damping, diameter]`．
`n_points` は節点数（= `n_segments + 1`）．これは BEM 側
[bem/core/BEM_inputfile_reader.hpp](../../bem/core/BEM_inputfile_reader.hpp) の
`mooring_*` パーサと一致するので，BEM 係留入力をそのまま `cable_solver` でも
静的平衡確認用に解ける．

トップレベルの `gravity` 等の共通スカラーは単一ライン形式と同じ．
