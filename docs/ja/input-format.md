---
layout: default
title: 入力ファイル形式
lang: ja
---

[English](../input-format.html)

# 入力ファイル形式

すべての入力ファイルはJSON形式です。シミュレーションには`settings.json`と1つ以上の物体定義ファイルが必要です。

## ディレクトリ構成

```
my_case/
├── settings.json      # シミュレーション設定
├── water.json         # 自由表面の流体定義
├── tank.json          # 造波水槽の壁面
├── wavemaker.json     # 造波機境界（任意）
├── absorber.json      # 消波装置（任意）
└── float.json         # 浮体（任意）
```

## settings.json

```json
{
  "output_directory": "./output",
  "input_files": ["water.json", "tank.json", "wavemaker.json"],
  "dt": 0.05,
  "end_time_step": 1000,
  "ALE": "pseudo_quad",
  "ALE_period": 1
}
```

### 主要パラメータ

| キー | 型 | 説明 |
|------|------|------|
| `output_directory` | string | 出力ファイルのパス（入力ディレクトリからの相対パス） |
| `input_files` | array | 読み込む物体定義ファイルのリスト |
| `dt` | number | 時間刻み |
| `end_time_step` | integer | 最大時間ステップ数 |
| `ALE` | string | ALE手法：`"linear"`、`"pseudo_quad"` |
| `ALE_period` | integer | N ステップごとにALEを適用 |

## 物体定義ファイル

各物体ファイルは表面メッシュと境界条件を定義します。

```json
{
  "name": "water",
  "type": "Fluid",
  "objfile": ["../../obj/Goring1979/water.obj"]
}
```

### 必須キー

| キー | 型 | 説明 |
|------|------|------|
| `name` | string | 物体の一意な識別子 |
| `type` | string | 物体の種類（下表参照） |
| `objfile` | array | OBJメッシュファイルのパス |

### 物体の種類

| 種類 | 説明 |
|------|------|
| `Fluid` | 自由表面流体領域 |
| `RigidBody` | 6自由度剛体 |
| `SoftBody` | 変形体 |
| `Absorber` | 消波境界 |

### 任意キー

| キー | 型 | 説明 |
|------|------|------|
| `ignore` | bool | `true`でこの物体をスキップ |
| `velocity` | array | 初速度 `[vx, vy, vz]` |
| `angular_velocity` | array | 初期角速度 |
| `center_of_mass` | array | 重心位置 |
| `mass` | number | 質量 |
| `MOI` | array | 慣性モーメント |
