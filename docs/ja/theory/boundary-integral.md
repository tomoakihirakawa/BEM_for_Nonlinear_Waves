---
layout: default
title: 境界積分方程式
lang: ja
---

[English version](../../theory/boundary-integral.html)

# 境界積分方程式

本ページでは、ポテンシャル流れの支配方程式から境界積分方程式を導出し、線形三角要素による離散化、係数行列の構成、多重節点の処理について解説します。

## 支配方程式

非圧縮・非回転流体の速度場は速度ポテンシャル $\phi$ を用いて $\mathbf{u} = \nabla\phi$ と表されます。連続の式から、$\phi$ は流体領域 $\Omega$ 内でラプラス方程式を満たします。

$$\nabla^2 \phi = 0 \quad \text{in } \Omega$$

## グリーンの定理からBIEへ

グリーンの第二恒等式を自由空間グリーン関数

$$G(\mathbf{x}, \mathbf{y}) = \frac{1}{4\pi |\mathbf{x} - \mathbf{y}|}$$

に対して適用すると、境界積分方程式が得られます。

$$c(\mathbf{x})\,\phi(\mathbf{x}) = \int_{\Gamma} \left[ G(\mathbf{x}, \mathbf{y})\,\frac{\partial \phi}{\partial n}(\mathbf{y}) - \phi(\mathbf{y})\,\frac{\partial G}{\partial n_y}(\mathbf{x}, \mathbf{y}) \right] dS(\mathbf{y})$$

ここで $c(\mathbf{x})$ は照合点 $\mathbf{x}$ における立体角係数であり、$\mathbf{x}$ が滑らかな境界上にある場合は $c = 1/2$ となります。$\Gamma$ は流体領域の全境界面を表します。

## 線形三角要素による離散化

境界面 $\Gamma$ を $N_e$ 個の三角要素で離散化します。各三角要素上で $\phi$ と $\partial\phi/\partial n$ を形状関数で補間します。

![線形三角要素](../../images/bem/schematic_linear_triangle_element.png)

### 形状関数

線形三角要素では、面積座標 $(\xi_1, \xi_2, \xi_3)$（ $\xi_3 = 1 - \xi_1 - \xi_2$ ）を用いて、要素内の任意の点における物理量を3頂点の値から線形補間します。

$$\phi(\xi_1, \xi_2) = \sum_{k=1}^{3} N_k(\xi_1, \xi_2)\,\phi_k$$

ここで $N_1 = \xi_1$, $N_2 = \xi_2$, $N_3 = 1 - \xi_1 - \xi_2$ です。

### ヤコビアン

パラメータ空間 $(\xi_1, \xi_2)$ から物理空間への写像のヤコビアンは次式で与えられます。

$$J = \left| \frac{\partial \mathbf{x}}{\partial \xi_1} \times \frac{\partial \mathbf{x}}{\partial \xi_2} \right|$$

面積分の変換に用いられ、$dS = J\,d\xi_1\,d\xi_2$ となります。数値積分にはDunavant求積法を使用します。

## 係数行列の構成

BIEを全ての境界節点 $\mathbf{x}_i$ （ $i = 1, \ldots, N$ ）について離散化すると、次の行列方程式が得られます。

$$\mathbf{M}\,\boldsymbol{\phi} = \mathbf{N}\,\frac{\partial \boldsymbol{\phi}}{\partial n}$$

ここで：

- **M行列**: $\frac{\partial G}{\partial n}$ の積分から構成される係数行列（対角成分に $c(\mathbf{x}_i)$ を含む）
- **N行列**: $G$ の積分から構成される係数行列

各要素について、Dunavant数値積分によりM行列およびN行列の成分を計算します。照合点を含む要素（特異積分）に対しては、特異点除去の手法を適用します。

### リジッドモードテクニック

M行列の対角成分 $M_{ii}$ の直接計算は特異積分を伴い精度確保が困難です。リジッドモードテクニックでは、 $\phi = 1$（一様ポテンシャル）を代入すると $\partial\phi/\partial n = 0$ となることを利用し、

$$M_{ii} = -\sum_{j \neq i} M_{ij}$$

として間接的に対角成分を求めます。これにより特異積分の直接計算を回避できます。

## GMRESとFMMの関連

境界値問題の解法として、直接法（LU分解）ではなく反復法のGMRES（Generalized Minimal Residual Method）を採用しています。

GMRESの各反復で必要な行列-ベクトル積の評価に[FMM](fmm.html)を適用することで、計算量を $O(N^2)$ から $O(N)$ に削減します。これにより、係数行列を陽的に保持する必要がなくなり、メモリ使用量も大幅に削減されます。

## 境界タイプ

![境界タイプの模式図](../../images/bem/schematic_boundary_types_without_float.png)

境界面の各節点には、以下の境界タイプが設定されます。

| 境界タイプ | 既知量 | 未知量 | 対応する境界 |
|-----------|--------|--------|-------------|
| Dirichlet | $\phi$ | $\partial\phi/\partial n$ | 自由表面 |
| Neumann | $\partial\phi/\partial n$ | $\phi$ | 固体壁面、海底 |

### 境界タイプの決定方法

各節点の境界タイプは、その節点が属する面の境界条件に基づいて決定されます。

- 全ての隣接面がディリクレ面の場合: **Dirichlet**節点
- 全ての隣接面がノイマン面の場合: **Neumann**節点
- ディリクレ面とノイマン面の両方に隣接する場合: **CORNER**節点

## 多重節点（CORNER節点）

CORNER節点は、ディリクレ境界とノイマン境界が交わる稜線上に位置する節点です。このような節点では $\partial\phi/\partial n$ が不連続となるため、物理的に1つの節点を複数の計算節点として扱います。

具体的には、CORNER節点ではノイマン面ごとに異なる $\partial\phi/\partial n$ の値を持たせます。これにより、1つの物理節点が複数の未知数を持つことになります。BIEの離散化においては、各計算節点に対して別々の方程式を立てます。

CORNER節点の処理は、波と構造物の相互作用問題において重要です。例えば、自由表面と固体壁面が交わる位置や、造波機表面と自由表面の接触線がこれに該当します。

## 関連ページ

- [理論概要](index.html)
- [浮体動力学](floating-body.html) --- 浮体がある場合の境界条件
- [高速多重極展開法](fmm.html) --- BIE評価の高速化
