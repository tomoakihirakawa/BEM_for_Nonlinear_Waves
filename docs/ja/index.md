---
layout: default
title: ホーム
---

[🇬🇧 English](../)

# 非線形波浪の境界要素法ソルバー

非線形自由表面波問題のための境界要素法（BEM）ソルバーです。

## 概要

本ソフトウェアは、境界要素法を用いて非線形自由表面境界条件を持つポテンシャル流れ問題を解きます。時間領域シミュレーションにはMixed Eulerian-Lagrangian（MEL）法を実装し、周波数領域解析にも対応しています。

### 主な機能

- **非線形造波** — ピストン/フラップ式造波機、孤立波、不規則波
- **波浪と物体の相互作用** — 6自由度剛体運動を伴う浮体
- **高速多重極展開法（FMM）** — O(N)の境界積分評価
- **ALEメッシュ管理** — ラプラシアン平滑化による適応的リメッシュ
- **係留索の動力学** — 集中質量ケーブルモデル
- **Metal GPU高速化** — Apple SiliconでのM2L変換

## はじめに

ビルド手順は[はじめに](getting-started.html)ガイドをご覧ください。

## ドキュメント

- [はじめに](getting-started.html) — ビルドと最初のシミュレーション実行
- [理論](theory.html) — 数学的定式化
- [入力ファイル形式](input-format.html) — JSON入力ファイルリファレンス

## 計算例

- [Goring (1979)](examples/goring1979.html) — 孤立波の生成と伝播

## ライセンス

LGPL-3.0-or-later。[LICENSE](https://github.com/tomoakihirakawa/BEM_for_Nonlinear_Waves/blob/main/LICENSE)を参照。
