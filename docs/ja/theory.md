---
layout: default
title: 理論
---

[🇬🇧 English](../theory.html)

# 理論的背景

## ポテンシャル流れの定式化

流体は非粘性、非圧縮、非回転と仮定します。速度場は速度ポテンシャル $\phi$ で記述され、ラプラス方程式を満たします：

$$\nabla^2 \phi = 0$$

## 境界積分方程式

グリーンの第二恒等式を用いると、境界積分方程式は次のようになります：

$$c(\mathbf{x})\phi(\mathbf{x}) = \int_S \left[ G(\mathbf{x}, \mathbf{y}) \frac{\partial \phi}{\partial n}(\mathbf{y}) - \phi(\mathbf{y}) \frac{\partial G}{\partial n}(\mathbf{x}, \mathbf{y}) \right] dS(\mathbf{y})$$

ここで $G(\mathbf{x}, \mathbf{y}) = \frac{1}{4\pi|\mathbf{x} - \mathbf{y}|}$ は自由空間グリーン関数、$c(\mathbf{x})$ は立体角係数です。

## Mixed Eulerian-Lagrangian（MEL）法

自由表面はラグランジュ的枠組みで追跡します。各時間ステップで：

1. **境界値問題の求解**: 自由表面上の $\phi$ と固体境界上の $\partial\phi/\partial n$ が与えられた条件で、境界積分方程式を解いて未知量を求めます。
2. **速度の計算**: 自由表面上の $\nabla\phi$ を計算します。
3. **自由表面の更新**: 運動学的・力学的境界条件を用いて自由表面の位置とポテンシャルを時間発展させます：

$$\frac{D\mathbf{x}}{Dt} = \nabla\phi, \quad \frac{D\phi}{Dt} = -gz + \frac{1}{2}|\nabla\phi|^2$$

## 高速多重極展開法（FMM）

FMMは境界積分の評価をO(N²)からO(N)に高速化します。本実装では以下を使用します：

- 球面調和関数展開による多重極係数・局所係数
- Multipole-to-Local（M2L）変換演算子
- 適応的八分木空間分割

### 利用可能なM2L手法

| 手法 | 説明 |
|------|------|
| `SimpleM2L` | 球面調和関数による直接変換（デフォルト） |
| `FourierM2L_double` | フーリエ基底回転、倍精度 |
| `FourierM2L_DoubleDouble` | フーリエ基底回転、倍々精度 |
| `PlaneWaveM2L` | 平面波展開法 |

## 任意ラグランジュ-オイラー（ALE）メッシュ移動

メッシュの歪みを防ぐため、ALE手法を用います：

1. ラグランジュ節点速度を計算
2. ラプラシアン平滑化で節点を再配置
3. ALE速度に対する物質微分を補正

## 造波

対応する造波手法：

- **ピストン/フラップ式造波機** — 規定運動
- **孤立波** 生成（Goring, 1979）
- **規則波** — 線形/2次造波理論
- **吸収境界** — 数値ビーチ
