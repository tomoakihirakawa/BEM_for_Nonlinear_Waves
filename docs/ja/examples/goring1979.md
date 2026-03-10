---
layout: default
title: "計算例: Goring (1979)"
---

[🇬🇧 English](../../examples/goring1979.html)

# 計算例：孤立波の生成（Goring, 1979）

ピストン式造波機を用いた孤立波生成のデモです。Goring (1979) の手法に従っています。

## 参考文献

Goring, D.G. (1979). *Tsunamis — The propagation of long waves onto a shelf*. W.M. Keck Laboratory of Hydraulics and Water Resources, California Institute of Technology, Report No. KH-R-38.

## 設定

2次元数値造波水槽：
- 水槽長さ：約10 m
- 水深：0.3 m
- 左端にピストン式造波機
- 右端に消波装置

造波機の変位は、指定した波高の孤立波を生成するために孤立波理論に従います。

## 入力ファイル

入力ファイルは `bem/input_files/` 内のGoring1979ディレクトリにあります。

### 実行方法

```bash
cd bem/build
./main_time_domain ../../bem/input_files/Goring1979_DT0d03_.../
```

## 期待される結果

孤立波は水槽内を以下の特性で伝播するはずです：
- 理論値と一致する波高
- 最小限の後続振動
- 複数の波高計位置での実測値との良好な一致
