---
layout: default
title: "計算例: DeepCWind 浮体式プラットフォーム"
lang: ja
---

[English](../../examples/deepcwind.html)

# DeepCWind 浮体式プラットフォーム

OC4 DeepCWind セミサブ型浮体式プラットフォームの波浪-物体相互作用
シミュレーションの計算例です。浮体は非線形自由表面流れと連成した
6自由度剛体運動を行います。

## 参考文献

Robertson, A., Jonkman, J., Masciola, M., Song, H., Goupee, A., Coulling, A.,
and Luan, C. (2014). *Definition of the Semisubmersible Floating System for
Phase II of OC4*. National Renewable Energy Laboratory, Technical Report
NREL/TP-5000-60601.

## 概要

DeepCWind セミサブは洋上風力タービンを支持するために設計された3本柱の
浮体式プラットフォームです。以下の構成要素からなります。

- 三角形配置の3本のオフセットカラム
- 中央カラム（タワー基部）
- カラムを接続するポンツーンおよびクロスブレース
- 各オフセットカラムの基部に設置されたヒーブプレート

本BEMシミュレーションでは、プラットフォームを造波水槽内の単一剛体として
モデル化します。ソルバは以下を含む完全非線形波浪-物体相互作用を計算します。

- 非線形自由表面境界条件
- 6自由度剛体力学（サージ、スウェイ、ヒーブ、ロール、ピッチ、ヨー）
- 物体表面における静水圧および動水圧の積分
- 水槽境界での造波と消波

## 入力ファイルの構成

DeepCWind ケースの入力ファイルは `bem/input_files/` 内の `DeepCWind_*`
という名前のディレクトリに格納されています。典型的なケースには5つの入力
ファイルが含まれます。

### settings.json

```json
{
    "max_dt": 0.5,
    "end_time_step": 10000000,
    "end_time": 100,
    "element": "linear",
    "ALE": "pseudo_quad",
    "ALEPERIOD": "1",
    "output_directory": "./output",
    "input_files": [
        "tank.json",
        "water.json",
        "float.json",
        "absorber.json"
    ],
    "meshing_options": [
        "surface_flip"
    ]
}
```

| パラメータ | 説明 |
|-----------|------|
| `max_dt` | 最大時間刻み幅（秒）。浮体問題では0.5秒が一般的。 |
| `end_time` | シミュレーション終了時刻。定常振動に達するのに十分な時間を設定する。 |
| `ALE` | 浮体シミュレーションでは `pseudo_quad` を推奨。 |

### water.json

流体領域メッシュを定義します。メッシュは初期喫水線において浮体表面と
一致する必要があります。

```json
{
    "name": "water",
    "type": "Fluid",
    "objfile": "path/to/water.obj"
}
```

### tank.json

水槽壁面（底面および遠方境界）を剛体として定義します。

```json
{
    "name": "tank",
    "type": "RigidBody",
    "isFixed": true,
    "objfile": "path/to/tank.obj"
}
```

### float.json

質量特性と6自由度力学を持つ浮体を定義します。

```json
{
    "name": "float",
    "type": "RigidBody",
    "velocity": "floating",
    "mass": 13556760.999999998,
    "COM": [0, 0, 171.93],
    "MOI": [13947000000.0, 15552000000.0, 13692000000.0],
    "objfile": "path/to/float.obj"
}
```

| パラメータ | 説明 |
|-----------|------|
| `velocity` | `"floating"` で6自由度自由浮体力学を有効にする |
| `mass` | 浮体の総質量 (kg) |
| `COM` | 重心座標 [x, y, z] (メートル) |
| `MOI` | 慣性モーメント [Ixx, Iyy, Izz] (kg m^2) |
| `objfile` | 浮体の表面メッシュ |

質量、重心、慣性モーメントはプラットフォーム設計と整合している必要が
あります。メッシュファイルは没水面の形状を定義します。

### absorber.json

水槽壁面からの反射を防止する消波境界を定義します。

```json
{
    "name": "absorber",
    "type": "Absorber",
    "isFixed": true,
    "objfile": "path/to/absorber.obj",
    "wave_theory_L": [3.7, 150.0, 180.0, 0, 180.0]
}
```

| パラメータ | 説明 |
|-----------|------|
| `type` | `"Absorber"` で消波アルゴリズムを有効にする |
| `wave_theory_L` | 消波パラメータ: [波長, ...形状パラメータ] |

## 主要なシミュレーションパラメータ

DeepCWind シミュレーションを設定する際には、以下の点に注意が必要です。

- **時間刻み**: `max_dt` は 0.25--0.5 秒を使用する。小さい値は精度が
  向上するが計算時間が増加する。
- **ALEスキーム**: `pseudo_quad` は `linear` よりも大きな物体運動に対して
  良好なメッシュ品質を提供する。
- **メッシュ解像度**: 流体メッシュは喫水線付近で物体形状と波浪-物体相互
  作用を解像するのに十分な細かさが必要。
- **終了時刻**: 過渡応答が減衰し定常振動に達するのに十分な時間を確保する
  （通常50--100秒以上）。

## ビルドと実行

```bash
cd bem/build
cmake -DSOURCE_FILE=../main_time_domain.cpp ../..
make -j$(sysctl -n hw.logicalcpu)
./main ../../bem/input_files/DeepCWind_DT0d5_ELEMlinear_ALEpseudo_quad_ALEPERIOD1/
```

入力ディレクトリのパスは実行したい DeepCWind ケースに応じて変更してください。

## 期待される結果

- 浮体は静水圧平衡位置に落ち着く
- 波浪荷重を受けてプラットフォームは6自由度の振動運動を示す
- ヒーブ、ピッチ、サージが支配的な応答モードとなる
- 物体表面の流体力が出力から抽出可能
- 消波装置が水槽境界での放射波を効果的に減衰させる

## 関連項目

- [浮体理論](../../theory/floating-body.html)
