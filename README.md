# BEM for Nonlinear Waves

## 概要

`BEM for Nonlinear Waves` は，非線形自由表面波問題を対象とする境界要素法（BEM）ベースの研究用コードである．
時間領域では Mixed Eulerian-Lagrangian（MEL）法に基づき，造波，波浪伝播，浮体運動，接触判定，自由表面節点の移動・再配置を扱う．
計算幾何学寄りのリメッシュ処理と時間領域 BEM を統合し，変形する自由表面・物体境界の品質を保ちながら，比較的安定した計算を実現することを重視している．
周波数領域では，浮体まわりの線形応答，放射問題，波強制問題，QTF 関連の計算を扱うための実装を進めている．
係留索や橋梁ケーブルは，`LumpedCable` と `cable_solver` による集中質量モデルとして分離している．
入力作成や結果確認のための GUI も含み，CLI と GUI の両方から利用できる．

## 主な機能

- **三角要素ベースの時間領域非線形 BEM**．任意形状や鋭利エッジ（薄板，ヒーブプレート等）を，四辺形要素より柔軟に扱える．
- **曲率適応リメッシュと BEM の統合**．曲面忠実度 θ（辺長 × 曲率）に基づく `split` / `collapse` / `flip` / `smooth` を，「試行してスコア改善するなら採用」という Hoppe 拡張方式で実行する．変形する自由表面・物体境界で長時間の安定性を保つ．
- **2 次曲面を介した物理量 transfer**．リメッシュ時，φ と φ_t を `TriShape<6>` の 2 次曲面上で Newton 精密投影により補間する．高次性を維持したままメッシュを再構成する．
- **個別メッシュからの自動境界条件設定**．流体・構造物・浮体を別々の OBJ として投入し，自動接触判定で `Dirichlet` / `Neumann` / `CORNER` / `multipleNode` を割り当てる．ユーザはメッシュを用意するだけで計算を開始できる．
- **多体連成と接触保護**．浮体 6DOF，係留索（`LumpedCable`），浮体‐浮体・浮体‐壁の衝突検出と保護領域を統合的に扱う．
- **動的境界属性変更**．時間発展中に `CORNER` / `multipleNode` が生成・消滅するケースに対応する．
- **FMM + Metal acceleration（任意）**．大規模積分を高速化する．Apple Silicon では GPU 経由の M2L 評価が可能である．

## まず読むもの

| 目的 | 読むもの |
|---|---|
| ビルド方法を確認する | [docs/ja/build.md](docs/ja/build.md) |
| 最初に BEM を動かす | [docs/ja/getting-started.md](docs/ja/getting-started.md) |
| BEM 入力 JSON を確認する | [docs/ja/input-format.md](docs/ja/input-format.md) |
| BEM GUI を使う | [docs/ja/gui.md](docs/ja/gui.md) |
| Goring (1979) の例題を見る | [docs/ja/examples/goring1979.md](docs/ja/examples/goring1979.md) |
| BEM 実装メモを見る | [bem/README.md](bem/README.md) |
| cable solver を使う | [cable/README.md](cable/README.md) |
| cable 入力 JSON を確認する | [cable/docs/input_json_schema.md](cable/docs/input_json_schema.md) |
| cable GUI を使う | [cable/gui/README.md](cable/gui/README.md) |

## 使い方の入口

### CLI

| 用途 | 実行ファイル | 主な入力 |
|---|---|---|
| 時間領域 BEM | `main_time_domain` | BEM 入力ディレクトリ |
| 周波数領域 BEM | `main_freq_domain` | BEM 入力ディレクトリ，周波数，DOF |
| 係留・ケーブル | `cable_solver` | cable JSON，出力ディレクトリ |

### GUI

| 用途 | 場所 | 内容 |
|---|---|---|
| BEM 入力 GUI | [gui/](gui/) | 入力ファイルの作成，編集，3D 表示 |
| cable GUI | [cable/gui/](cable/gui/) | `cable_solver` の入力作成，実行，結果表示 |

## リポジトリ構成

```text
BEM_for_Nonlinear_Waves/
├── bem/
│   ├── core/              # BEM 共通処理，入力読み込み，境界条件，BVP
│   ├── time_domain/       # 時間領域 BEM ソルバ
│   ├── frequency_domain/  # 周波数領域 BEM ソルバ
│   ├── input_files/       # BEM 入力例
│   ├── metal_m2l/         # Metal M2L acceleration
│   └── metal_nearfield/   # Metal nearfield acceleration
├── lib/
│   ├── include/           # 共有ヘッダ，Network，FMM，幾何，ODE
│   └── src/               # 共有実装
├── cable/                 # 係留索・ケーブルソルバ，GUI，例題
├── gui/                   # BEM 入力 GUI
├── docs/                  # ドキュメント
├── obj/                   # メッシュ・ベンチマーク形状
├── third_party/tetgen/    # TetGen header + pre-built library
└── CMakeLists.txt         # 共通 CMake エントリ
```

## 最小ビルド

詳細は [docs/ja/build.md](docs/ja/build.md) を参照する．代表的な入口は次の 3 つである．

```bash
mkdir -p builds/build_bem_time
cd builds/build_bem_time
cmake -DCMAKE_BUILD_TYPE=Release -DSOURCE_FILE=bem/time_domain/main_time_domain.cpp ../..
cmake --build . -j$(sysctl -n hw.logicalcpu)
```

```bash
mkdir -p builds/build_bem_freq
cd builds/build_bem_freq
cmake -DCMAKE_BUILD_TYPE=Release -DSOURCE_FILE=bem/frequency_domain/main_freq_domain.cpp ../..
cmake --build . -j$(sysctl -n hw.logicalcpu)
```

```bash
mkdir -p cable/build_solver
cd cable/build_solver
cmake -DCMAKE_BUILD_TYPE=Release -DSOURCE_FILE=cable/cable_solver.cpp -DOUTPUT_NAME=cable_solver ../..
cmake --build . -j$(sysctl -n hw.logicalcpu)
```

## 最小実行

### CLI

```bash
./main_time_domain /path/to/bem/input_directory
```

```bash
./main_freq_domain /path/to/bem/input_directory --omega 0.5 0.8 1.1
```

```bash
./cable_solver cable/gui/examples/synthetic/catenary_500m.json /tmp/cable_out/
```

### GUI

```bash
cd gui
./run.sh
```

```bash
cd cable/gui
./run.sh
```

## モジュール別ドキュメント

| モジュール | ドキュメント |
|---|---|
| ビルド | [docs/ja/build.md](docs/ja/build.md) |
| BEM 全体 | [bem/README.md](bem/README.md) |
| BEM 入力形式 | [docs/ja/input-format.md](docs/ja/input-format.md) |
| BEM 例題 | [bem/EXAMPLES.md](bem/EXAMPLES.md) |
| Goring (1979) 例題 | [docs/ja/examples/goring1979.md](docs/ja/examples/goring1979.md) |
| 周波数領域 BEM | [bem/docs/周波数領域BEMの概要.md](bem/docs/周波数領域BEMの概要.md) |
| Vortex Particle Method | [bem/docs/VPMの概要.md](bem/docs/VPMの概要.md) |
| cable solver | [cable/README.md](cable/README.md) |
| cable 入力形式 | [cable/docs/input_json_schema.md](cable/docs/input_json_schema.md) |
| cable GUI | [cable/gui/README.md](cable/gui/README.md) |
| cable 例題 | [cable/gui/examples/README.md](cable/gui/examples/README.md) |
| BEM GUI | [gui/README.md](gui/README.md) |
| TetGen | [third_party/tetgen/README.md](third_party/tetgen/README.md) |

## 開発状況

| 項目 | 状態 |
|---|---|
| 時間領域 BEM | 主開発対象 |
| 周波数領域 BEM | 開発中 |
| cable solver | 独立 CLI として利用可能 |
| cable GUI | macOS / Apple Silicon 中心に検証 |
| BEM GUI | 入力作成・可視化用として開発中 |
| Metal acceleration | 任意・実験的 |

## ライセンス

このプロジェクトは [LGPL-3.0-or-later](LICENSE) で公開する．

一部の third-party component は個別のライセンスを持つ．

- TetGen: AGPL-3.0，詳細は [third_party/tetgen/LICENSE](third_party/tetgen/LICENSE) を参照．
- ankerl::unordered_dense: MIT．
- pdqsort: zlib．

TetGen は四面体分割機能のための任意依存である．TetGen を使わない場合は `-DUSE_TETGEN=OFF` を指定する．

## 引用

研究で利用する場合は，次を引用する．

```bibtex
@software{hirakawa2025bem,
  author = {Hirakawa, Tomoaki},
  title = {BEM for Nonlinear Waves},
  year = {2025},
  url = {https://github.com/tomoakihirakawa/BEM_for_Nonlinear_Waves}
}
```
