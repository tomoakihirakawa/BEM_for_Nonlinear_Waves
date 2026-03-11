---
layout: default
title: はじめに
lang: ja
---

[English](../getting-started.html)

# はじめに

## 必要環境

| 要件 | 最小バージョン | 備考 |
|------|---------------|------|
| C++コンパイラ | GCC 12+ または Clang 15+ | `__float128`サポートのためGCC推奨 |
| CMake | 3.16+ | |
| LAPACK/BLAS | — | 線形代数計算に必要 |
| OpenMP | — | 任意だが推奨 |

### macOS (Homebrew)

```bash
brew install gcc cmake lapack
```

### Ubuntu/Debian

```bash
sudo apt install g++ cmake liblapack-dev libblas-dev libgomp1
```

## ビルド

### リポジトリのクローン

```bash
git clone https://github.com/tomoakihirakawa/BEM_for_Nonlinear_Waves.git
cd BEM_for_Nonlinear_Waves
```

### 時間領域ソルバー

```bash
mkdir -p bem/build && cd bem/build
cmake -DSOURCE_FILE=../main_time_domain.cpp ../..
cmake --build . -j$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)
```

### 周波数領域ソルバー

```bash
mkdir -p build-freq && cd build-freq
cmake -DSOURCE_FILE=bem/main_freq_domain.cpp ..
cmake --build . -j$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)
```

## シミュレーションの実行

```bash
./main_time_domain /path/to/input_files/
```

入力ディレクトリには以下が必要です：
- `settings.json` — シミュレーションパラメータ（時間刻み、計算時間など）
- 1つ以上の物体/流体定義ファイル（例：`water.json`、`tank.json`）

各物体ファイルでは以下を指定します：
- `name` — 識別子
- `type` — `Fluid`、`RigidBody`、`SoftBody`など
- `objfile` — 表面メッシュファイルのパス（OBJ形式）

## CMakeオプション

```bash
cmake -DSOURCE_FILE=../main_time_domain.cpp \
      -DBEM_COMPILER=gcc \
      -DUSE_TETGEN=ON \
      -DFMM_M2L_METHOD=SimpleM2L \
      ../..
```

すべてのオプションは[README](https://github.com/tomoakihirakawa/BEM_for_Nonlinear_Waves#cmake-options)を参照してください。

## TetGenなしでのビルド

TetGenはAGPL-3.0ライセンスです。TetGenなしでビルドするには：

```bash
cmake -DUSE_TETGEN=OFF -DSOURCE_FILE=../main_time_domain.cpp ../..
```

四面体分割機能は無効になりますが、BEMソルバーは表面のみの計算で完全に動作します。
