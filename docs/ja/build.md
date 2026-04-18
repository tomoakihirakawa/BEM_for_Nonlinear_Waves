---
layout: default
title: ビルド方法
lang: ja
---

# ビルド方法

## 概要

このページでは，`BEM for Nonlinear Waves` の主要な実行ファイルをビルドする方法をまとめる．
対象は時間領域 BEM，周波数領域 BEM，`cable_solver`，GUI である．
ルートの `CMakeLists.txt` は `SOURCE_FILE` に指定したソースを入口として実行ファイルを作る．
コンパイラや CMake option を変える場合，CMake cache の影響を避けるため，新しいビルドディレクトリを使う．

## 必要環境

| 項目 | 目安 | 備考 |
|---|---|---|
| C++ compiler | GCC 12+ または Clang 15+ | C++23 を使う．GCC は `__float128` と `quadmath` を使える |
| CMake | 3.16+ | ルート `CMakeLists.txt` を使う |
| LAPACK / BLAS | 必須 | 線形代数計算に使う |
| OpenMP | 推奨 | 並列化に使う |
| Python | GUI 用 | `gui/` と `cable/gui/` で使う |

macOS + Homebrew の例:

```bash
brew install gcc cmake lapack openblas libomp
```

Ubuntu / Debian の例:

```bash
sudo apt install g++ cmake liblapack-dev libblas-dev libgomp1
```

## ビルドディレクトリの方針

ビルド成果物は，ソースツリー直下の実装ディレクトリではなく，専用のビルドディレクトリに分ける．
推奨する基本形は次である．

| 対象 | 推奨ビルドディレクトリ |
|---|---|
| 時間領域 BEM | `builds/build_bem_time/` |
| 周波数領域 BEM | `builds/build_bem_freq/` |
| cable solver | `cable/build_solver/` |

`SOURCE_FILE` はルート `CMakeLists.txt` から見た相対パスで指定する．
`OUTPUT_NAME` を指定しない場合，実行ファイル名は `SOURCE_FILE` の basename から決まる．

## 時間領域 BEM

```bash
mkdir -p builds/build_bem_time
cd builds/build_bem_time
cmake -DCMAKE_BUILD_TYPE=Release -DSOURCE_FILE=bem/time_domain/main_time_domain.cpp ../..
cmake --build . -j$(sysctl -n hw.logicalcpu)
```

生成される実行ファイル:

```text
builds/build_bem_time/main_time_domain
```

実行例:

```bash
./main_time_domain /path/to/bem/input_directory
```

入力ディレクトリには `settings.json` と，`settings.json` から参照される物体・流体定義 JSON が必要である．

## 周波数領域 BEM

```bash
mkdir -p builds/build_bem_freq
cd builds/build_bem_freq
cmake -DCMAKE_BUILD_TYPE=Release -DSOURCE_FILE=bem/frequency_domain/main_freq_domain.cpp ../..
cmake --build . -j$(sysctl -n hw.logicalcpu)
```

生成される実行ファイル:

```text
builds/build_bem_freq/main_freq_domain
```

実行例:

```bash
./main_freq_domain /path/to/bem/input_directory --omega 0.5 0.8 1.1
```

DOF を指定する場合:

```bash
./main_freq_domain /path/to/bem/input_directory --omega 0.5 0.8 1.1 --dofs 0 1 2
```

## cable solver

```bash
mkdir -p cable/build_solver
cd cable/build_solver
cmake -DCMAKE_BUILD_TYPE=Release -DSOURCE_FILE=cable/cable_solver.cpp -DOUTPUT_NAME=cable_solver ../..
cmake --build . -j$(sysctl -n hw.logicalcpu)
```

生成される実行ファイル:

```text
cable/build_solver/cable_solver
```

実行例:

```bash
./cable_solver ../gui/examples/synthetic/catenary_500m.json /tmp/cable_out/
```

`cable_solver` は単一ケーブル，複数ケーブル，BEM 互換の係留入力を自動判別する．
入力 JSON の詳細は [cable/docs/input_json_schema.md](../../cable/docs/input_json_schema.md) を参照する．

## GUI

### BEM GUI

BEM 入力 GUI は [gui/](../../gui/) にある．
依存関係と既知の互換性は [gui/README.md](../../gui/README.md) を参照する．

```bash
cd gui
./run.sh
```

### cable GUI

cable GUI は [cable/gui/](../../cable/gui/) にある．
内部では `cable_solver` を呼び出すため，先に `cable/build_solver/cable_solver` をビルドしておく．

```bash
cd cable/gui
./run.sh
```

詳細は [cable/gui/README.md](../../cable/gui/README.md) を参照する．

## CMake オプション

| Option | Default | 説明 |
|---|---|---|
| `BEM_COMPILER` | `gcc` | `gcc` または `clang` を選ぶ |
| `BEM_ENABLE_OPENMP` | `ON` | OpenMP を有効化する |
| `BEM_DISABLE_FLOAT128` | `OFF` | `std::float128_t` / `quadmath` 系のコードパスを無効化する |
| `FMM_M2L_METHOD` | `SimpleM2L` | FMM M2L translation method を選ぶ |
| `USE_TETGEN` | `ON` | TetGen を使う |
| `USE_METAL_M2L` | `OFF` | Metal M2L acceleration を使う |
| `USE_METAL_NEARFIELD` | `OFF` | Metal nearfield acceleration を使う |
| `CCACHE_ENABLE` | `ON` | `ccache` があれば使う |

Clang を使う例:

```bash
cmake -DCMAKE_BUILD_TYPE=Release \
      -DBEM_COMPILER=clang \
      -DSOURCE_FILE=bem/time_domain/main_time_domain.cpp \
      ../..
```

CMake は最初に選んだ compiler を cache する．compiler を切り替える場合は，新しいビルドディレクトリを使う．

## TetGen なしでのビルド

TetGen は AGPL-3.0 licensed component であり，四面体分割機能のために任意で使う．
TetGen を使わずにビルドする場合:

```bash
cmake -DCMAKE_BUILD_TYPE=Release \
      -DUSE_TETGEN=OFF \
      -DSOURCE_FILE=bem/time_domain/main_time_domain.cpp \
      ../..
```

TetGen の同梱物とライセンスは [third_party/tetgen/README.md](../../third_party/tetgen/README.md) と [third_party/tetgen/LICENSE](../../third_party/tetgen/LICENSE) を参照する．

## Metal acceleration

Metal acceleration は Apple Silicon 向けの任意・実験的な機能である．
`USE_METAL_M2L` または `USE_METAL_NEARFIELD` を有効にする場合，対応する Metal library を先にビルドしておく必要がある．

M2L:

```bash
cd bem/metal_m2l
mkdir -p build
cd build
cmake ..
cmake --build .
```

その後，BEM 本体を次のように構成する．

```bash
cmake -DCMAKE_BUILD_TYPE=Release \
      -DUSE_METAL_M2L=ON \
      -DSOURCE_FILE=bem/time_domain/main_time_domain.cpp \
      ../..
```

Nearfield:

```bash
cd bem/metal_nearfield
mkdir -p build
cd build
cmake ..
cmake --build .
```

その後，BEM 本体を次のように構成する．

```bash
cmake -DCMAKE_BUILD_TYPE=Release \
      -DUSE_METAL_NEARFIELD=ON \
      -DSOURCE_FILE=bem/time_domain/main_time_domain.cpp \
      ../..
```

## トラブルシューティング

### compiler を切り替えても反映されない

CMake cache が残っている可能性が高い．新しいビルドディレクトリを作る．

### `lapack` または `blas` が見つからない

macOS では Homebrew の `lapack`，`openblas` を入れる．
必要なら `LDFLAGS` と `CPPFLAGS` で Homebrew の install path を明示する．

### Clang で OpenMP が見つからない

Homebrew の `libomp` または LLVM を入れる．

```bash
brew install llvm libomp
```

### `cable_solver` が GUI から見つからない

`cable/build_solver/cable_solver` をビルドする．
別の場所に置く場合は，`PYCABLE_SOLVER_PATH` で実行ファイルを指定する．

```bash
export PYCABLE_SOLVER_PATH=/path/to/cable_solver
```
