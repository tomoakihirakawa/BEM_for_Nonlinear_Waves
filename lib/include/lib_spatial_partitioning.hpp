#pragma once

#include <algorithm>
#include <execution>
#include <limits>
#include <ranges>

#include "lib_Fourier.hpp"
#include "lib_multipole_expansion.hpp"
#include "moments.hpp" // 外部化した Moments

namespace spatial_partitioning_detail {
template <class It, class Pred>
inline bool any_of_unseq(It first, It last, Pred pred) {
#if defined(__cpp_lib_execution)
  return std::any_of(std::execution::unseq, first, last, pred);
#else
  return std::any_of(first, last, pred);
#endif
}

template <class It, class Pred>
inline bool all_of_unseq(It first, It last, Pred pred) {
#if defined(__cpp_lib_execution)
  return std::all_of(std::execution::unseq, first, last, pred);
#else
  return std::all_of(first, last, pred);
#endif
}

template <class It, class Fn>
inline void for_each_unseq(It first, It last, Fn fn) {
#if defined(__cpp_lib_execution)
  std::for_each(std::execution::unseq, first, last, fn);
#else
  std::for_each(first, last, fn);
#endif
}

template <class It, class Fn>
inline void for_each_par_unseq(It first, It last, Fn fn) {
#if defined(__cpp_lib_execution)
  std::for_each(std::execution::par_unseq, first, last, fn);
#else
  std::for_each(first, last, fn);
#endif
}
} // namespace spatial_partitioning_detail

/*DOC_EXTRACT lib_spatial_partitioning

## 空間分割（space_partitioning） `Buckets<T, N>`

## 概要

`Buckets<T, N>` は、3次元空間内のオブジェクトを効率的に管理、検索するための汎用的な空間分割データ構造です。主な目的は、**高速な近傍探索**と、**高速多重極法 (FMM: Fast Multipole Method)** のための階層的なデータ管理基盤を提供することです。

このクラスは、空間を均一な格子（セル/バケット）に分割するだけでなく、オブジェクトの密度に応じて各セルを再帰的に細分化し、オクツリーのような**階層ツリー構造**を動的に構築できます。

**主な機能:**

  * **フラットな空間分割**: 特定の座標周辺のオブジェクトを高速にリストアップします。
  * **階層的なツリー構築**: オブジェクトの分布に応じて空間を再帰的に分割し、計算の局所性を高めます。
  * **FMMサポート**: 多重極展開と局所展開のモーメント、および相互作用リストを各ノードで保持し、FMM計算をサポートします。
  * **パフォーマンス最適化**: 更新性能と走査性能を両立させるためのデータ構造や、C++17の並列アルゴリズムを活用した高速な操作を提供します。

### テンプレートパラメータ

| パラメータ | 説明 |
| :--- | :--- |
| `T` | 格納するオブジェクトの型。通常はオブジェクトへのポインタ（例: `Particle*`, `Node*`）を指定します。`T` 型のオブジェクト `p` は、`p->X` という形で `Tddd` 型の位置座標を公開している必要があります。 |
| `N` | FMMで使用するモーメントの最大次数。`Moments<T, N>` クラスの実装に依存します。`N=0` の場合はFMM関連の機能は実質的に無効化されます。 |

-----

## コアコンセプト

`Buckets` クラスは、2つの主要なモードで動作します。「フラットな空間分割」と「階層的なツリー構造」です。

### 1\. フラットな空間分割 (Flat Partition)

最も基本的な利用形態です。コンストラクタで指定された3次元のバウンディングボックスを、指定されたセル幅 `dL` で均一な格子に分割します。

  * **オブジェクトの追加**: 各オブジェクトは、その位置座標に基づいて対応するセル（バケット）に格納されます。
  * **近傍探索**: `apply(x, d, ...)` などの関数を呼び出すことで、座標 `x` から距離 `d` の範囲内にあるセル群を効率的に特定し、その中のオブジェクトに対して一括で処理を適用できます。これにより、全オブジェクトを走査するよりも劇的に計算コストを削減できます。

### 2\. 階層的なツリー構造 (Hierarchy)

`generateTree()` メソッドを呼び出すことで、フラットな構造から階層的なツリー構造へ移行します。

  * **再帰的分割**: 各セルは、内部のオブジェクト数が特定の条件（ラムダ式で指定可能）を満たす場合に、さらに小さな8つの子セル（3Dの場合）に再帰的に分割されます。これにより、オブジェクトが密集している領域は細かく、疎な領域は粗いままの、適応的な空間分割が実現されます。
  * **FMMの基盤**: このツリー構造はFMMの根幹をなします。各ノード（セル）が多重極モーメントや局所展開モーメントを保持し、親子関係や近傍関係を利用してM2M（Multipole-to-Multipole）、M2L（Multipole-to-Local）、L2L（Local-to-Local）といった変換を効率的に行います。

-----

## 内部データ構造

`Buckets` クラスは、パフォーマンスと柔軟性を両立させるために、複数のデータコンテナを内部で管理しています。

| メンバ変数 | 型 | 説明 |
| :--- | :--- | :--- |
| `data` | `vector<vector<vector<unordered_set<T>>>>` | オブジェクトを格納する主要な3Dコンテナ。**`unordered_set`** を使用しているため、オブジェクトの**追加・削除が平均 O(1)** と非常に高速です。 |
| `data_vector` | `vector<vector<vector<vector<T>>>>` | `data` の内容をコピーしたキャッシュ用の3Dコンテナ。**`vector`** はデータがメモリ上で連続しているため、範囲を指定した**走査（イテレーション）が `unordered_set` よりも高速**です。`setVector()` を呼び出すことで `data` と同期されます。 |
| `data1D` | `unordered_set<T>` | このバケット（およびすべての子孫バケット）に含まれる全オブジェクトの一意な集合。 |
| `children` | `vector<vector<vector<Buckets*>>>` | 階層構造における子ノード（子バケット）へのrawポインタを保持します。メモリ管理は手動で行われます。 |
| `level_buckets` | `vector<vector<Buckets*>>` | **ルートノードのみが保持**。ツリーの各レベルに存在する非空ノードへのポインタを格納した配列。`forEachAtLevel` などで特定のレベルのノードを高速に走査するために使われます。 |
| `deepest_level_buckets` | `vector<Buckets*>` | **ルートノードのみが保持**。ツリーの末端ノード（葉ノード）の集合。 |
| `Moments...` | `Moments<T,N>` | FMM用の多重極展開モーメントと局所展開モーメントを保持するオブジェクト。 |
| `buckets_for_M2L`, `buckets_near` | `vector<Buckets*>` | FMMの計算で必要となる相互作用リスト。 |

**不変条件**:

  * `vector_is_set == true` のとき、`data` と `data_vector` の内容は同一です。
  * 子ノードは、親ノードの領域を重複なく正確に分割します。

-----

## APIガイド

### A. 初期化と設定

  * `Buckets(bounds, dL)`: 指定された境界 `bounds` とセル幅 `dL` でルートバケットを初期化します。
  * `initialize(bounds, dL)`: 既存のインスタンスを再初期化します。
  * `setVector()`: `data` (`unordered_set`) の内容を `data_vector` (`vector`) にコピーします。大量の検索や走査処理の前に呼び出すことで、パフォーマンスが向上します。

### B. オブジェクトの管理 (追加・削除)

  * `add(x, p)`: 座標 `x` にオブジェクト `p` を追加します。ツリーが存在する場合は、対応する子孫ノードにも再帰的に追加されます。
  * `erase(p)`: オブジェクト `p` を削除します。ツリー全体から再帰的に削除されます。
  * `add_grow(x, p)`: オブジェクト `p` を追加し、分割条件を満たせばツリーを動的に成長（分割）させます。
  * `erase_shrink(p)`: オブジェクト `p` を削除し、その結果セルが空になった場合にノードを削除（縮小）します。
  * `clear()`: すべてのオブジェクトを削除します。

### C. 階層ツリーの操作

  * `generateTree(condition)`: 指定された `condition` (ラムダ式) を満たすノードを再帰的に分割し、階層ツリーを構築します。
      * 例: `B.generateTree([](auto b){ return b->data1D.size() > 8; });` // オブジェクト数が8を超えるバケットを分割
  * `forEachAll(func)`: ツリーの全ノード (ルート含む) に対して深さ優先で `func` を適用します。
  * `forEachAtLevel(level, func)`: 指定した `level` のすべてのノードに対して `func` を適用します。 (ルートが事前計算したリストを使用するため高速)
  * `forEachAtDeepest(func)`: ツリーのすべての葉ノードに対して `func` を適用します。
  * `getBucketAtLevel(level, x)`: 座標 `x` を含む、指定 `level` のノードへのポインタを返します。

### D. 近傍探索と走査

これらの関数は、内部で `vector_is_set` フラグをチェックし、`data_vector` が利用可能であればそちらを使って高速に処理を実行します。

  * `apply(x, d, func)`: 座標 `x` からの距離 `d` の範囲（立方体近似）に含まれる**すべてのオブジェクト**に `func` を適用します。
  * `any_of(x, d, pred)`: 範囲内に `pred` を満たすオブジェクトが**一つでも存在すれば** `true` を返します。
  * `all_of(x, d, pred)`: 範囲内の**すべてのオブジェクトが** `pred` を満たせば `true` を返します。
  * `none_of(x, d, pred)`: 範囲内に `pred` を満たすオブジェクトが**一つも存在しなければ** `true` を返します。

-----

## 代表的なワークフロー

```cpp
#include "lib_spatial_partitioning.hpp"

// オブジェクトの型 (例: networkPoint*)
using MyObject = networkPoint*;

// 1. ルートバケットの生成
CoordinateBounds bounds = ...; // 対象となる空間全体の境界
double dL = 1.0;               // 分割するセルの基本サイズ
Buckets<MyObject, 8> B(bounds, dL); // FMM次数8で初期化

// 2. オブジェクトの追加
for (MyObject p : all_points) {
    B.add(p->X, p);
}

// 3. (FMM利用時) 階層ツリーの生成
//    オブジェクト数が16を超えたノードを再帰的に分割する
B.generateTree([](auto b){
    return b->data1D.size() > 16;
});

// 4. (パフォーマンス向上) 走査用キャッシュの作成
//    これから多くの近傍探索を行う前に呼び出す
B.setVector();

// 5. 近傍クエリの実行
Tddd search_pos = {1.0, 2.0, 3.0};
double radius = 5.0;
B.apply(search_pos, radius, [](MyObject p){
    // 見つかったオブジェクト p に対する処理
    process(p);
});

// 6. (動的シミュレーション時) オブジェクトの移動とツリーの再構築
//    (batchRelocateAndMaybeRegrow のような外部ユーティリティを想定)
relocate_points();
B.rebuild_tree_if_needed();
```

-----

## 設計思想とトレードオフ

  * **更新性能 vs 走査性能**: `unordered_set` (`data`) と `vector` (`data_vector`) の二重持ちは設計の核です。前者はオブジェクトの頻繁な追加・削除に、後者は静的なオブジェクト群に対する大量の近傍クエリに強いという特性を両立させています。`setVector()` がその切り替えスイッチの役割を果たします。
  * **メモリ管理**: 子ノードを `raw pointer` (`Buckets*`) で管理することで、スマートポインタのオーバーヘッドを避け、パフォーマンスを優先しています。ただし、これによりメモリリークのリスクが生じるため、デストラクタや `erase_shrink` などで手動のメモリ管理が必須となります。将来的には `std::unique_ptr` への移行が検討されています。
  * **探索範囲の近似**: 球形の探索範囲は、実装の簡略化とSIMD命令の親和性を高めるため、軸並行境界ボックス（AABB）で近似されます。つまり、検索は少し大きめの立方体領域に対して行われます。厳密な球形範囲が必要な場合は、`apply` 内のラムダ式で追加の距離二乗チェックを行う必要があります。
  * **ルートノードへの情報集約**: `level_buckets` などの階層情報はルートノードに集約されます。これにより、特定のレベルのノード全体にアクセスする操作が O(対象ノード数) となり、ツリーを毎回トラバースするコストを削減しています。

## スレッド安全性

本クラスは、**スレッドセーフではありません**。

  * **書き込み操作** (`add`, `erase`, `generateTree` など) と、**読み取り操作** (`apply`, `any_of` など) を異なるスレッドで同時に実行すると、データ競合が発生し、未定義の動作を引き起こします。

安全な利用のためには、プログラムのフェーズを明確に分離する必要があります。

1.  **更新フェーズ**: シングルスレッドでオブジェクトの追加、削除、ツリーの再構築を行います。
2.  `setVector()` を呼び出します。
3.  **読み取りフェーズ**: `apply` などの読み取り専用操作を複数のスレッドで並列に実行します。内部で C++17 Parallel Algorithms や OpenMP が使われているため、読み取り操作自体は並列化の恩恵を受けます。

*/

// Buckets is derived from CoordinateBounds
// N <= 8, 10だと破綻することがあった．
//! ./fast ./input_files/Goring1979_DT0d05_MESHwater0d04refined.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
template <typename T /*Stored Object*/, int N = 0>
struct Buckets : public CoordinateBounds {

  static constexpr int expansion_order = N; // FMM展開次数を公開
  using sizeType = int;
  using ST = sizeType;
  using ST2 = std::array<sizeType, 2>;
  using ST3 = std::array<sizeType, 3>;
  using ST6 = std::array<sizeType, 6>;
  ST xsize, ysize, zsize;
  Tddd center;
  ST3 dn;
  std::vector<std::vector<std::vector<std::unordered_set<T>>>> data;
  std::vector<std::vector<std::vector<std::vector<T>>>> data_vector;
  bool vector_is_set;
  std::unordered_set<T> data1D;
  std::vector<T> data1D_vector;
  std::unordered_map<T, ST3> map_to_ijk;
  double dL;
  double inv_dL; // Inverse of dL for fast multiplication instead of division
  // Cached lower bounds for fast index calculation
  double xmin, ymin, zmin;
  double bucketVolume() const { return std::pow(this->dL, 3.); };

  Buckets() = default;
  Buckets(const CoordinateBounds& c_bounds) : CoordinateBounds(c_bounds) {
    auto [xmin, xmax] = c_bounds.bounds[0];
    auto [ymin, ymax] = c_bounds.bounds[1];
    auto [zmin, zmax] = c_bounds.bounds[2];
    double max_dL = std::max(std::abs(xmax - xmin), std::max(std::abs(ymax - ymin), std::abs(zmax - zmin)));
    initialize(this->bounds, max_dL);
    this->center = this->X;
  };
  Buckets(const CoordinateBounds& c_bounds, const double dL_IN) : CoordinateBounds(c_bounds) {
    initialize(this->bounds, dL_IN);
    this->center = this->X;
  };
  Buckets(const T3Tdd& boundingboxIN, const double dL_IN) : CoordinateBounds(boundingboxIN) {
    initialize(this->bounds, dL_IN);
    this->center = this->X;
  };

  ~Buckets() {
    this->data.clear();
    this->data_vector.clear();
    this->data1D.clear();
    this->data1D_vector.clear();
    this->map_to_ijk.clear();
    for (auto& i : this->children)
      for (auto& j : i)
        for (auto& k : j)
          if (k != nullptr) {
            delete k;
            k = nullptr;
          }
    this->children.clear();
  }

  std::vector<std::array<int, 3>> vector_ijk;

  ST calculate_size(const auto& bounds, const double& dL) {
    const std::size_t max_size = 1000000;
    double size = std::ceil((std::get<1>(bounds) - std::get<0>(bounds)) / dL);
    if (size <= 0.)
      size = 1.;
    if (size > max_size) {
      std::cerr << "The size of the bucket is too large (" << size << "). Please reduce the bucket size." << std::endl;
      exit(1);
    }
    return static_cast<ST>(size);
  };

  bool is_data_assigned_as_3D = false;
  void initialize(const auto& boundingboxIN, const double dL_IN) {

    this->vector_is_set = false;
    this->data1D.clear();
    this->data1D_vector.clear();
    this->map_to_ijk.clear();

    //

    CoordinateBounds::setBounds(boundingboxIN);
    this->dL = dL_IN;
    this->inv_dL = 1.0 / dL_IN; // Precompute inverse for fast multiplication
    // Cache lower bounds for fast index calculation
    this->xmin = std::get<0>(this->xbounds());
    this->ymin = std::get<0>(this->ybounds());
    this->zmin = std::get<0>(this->zbounds());

    auto xsize_last = this->xsize;
    auto ysize_last = this->ysize;
    auto zsize_last = this->zsize;

    this->xsize = calculate_size(this->xbounds(), this->dL);
    this->ysize = calculate_size(this->ybounds(), this->dL);
    this->zsize = calculate_size(this->zbounds(), this->dL);

    this->vector_is_set = false;

    if (this->is_data_assigned_as_3D && xsize == xsize_last && ysize == ysize_last && zsize == zsize_last) {
#pragma omp parallel for
      for (auto i = 0; i < this->xsize; ++i)
        for (auto j = 0; j < this->ysize; ++j)
          for (auto k = 0; k < this->zsize; ++k) {
            this->data[i][j][k].clear();
            this->data_vector[i][j][k].clear();
          }
      std::cout << Red << "successfully re-initialized Buckets with the same size (" << this->xsize << " x " << this->ysize << " x " << this->zsize << ")." << colorReset << std::endl;
    } else {
      vector_ijk.clear();
      vector_ijk.reserve(this->xsize * this->ysize * this->zsize);
      for (auto i = 0; i < this->xsize; ++i)
        for (auto j = 0; j < this->ysize; ++j)
          for (auto k = 0; k < this->zsize; ++k)
            vector_ijk.emplace_back(std::array<int, 3>{i, j, k});

      this->dn = {xsize, ysize, zsize};
      this->data.assign(xsize, std::vector<std::vector<std::unordered_set<T>>>(ysize, std::vector<std::unordered_set<T>>(zsize)));
      this->data_vector.assign(xsize, std::vector<std::vector<std::vector<T>>>(ysize, std::vector<std::vector<T>>(zsize)));
      this->is_data_assigned_as_3D = true;
    }

    this->children.reserve(this->xsize);
    /* -------------------------------------------------------------------------- */

    this->MomentsMultipoleExpansion.initialize(this->X);
    this->MomentsLocalExpansion.initialize(this->X);
  }

  //@ -------------------------------------------------------------------------- */
  //@                                  インデックス変換                             */
  //@ -------------------------------------------------------------------------- */

  /*
  std::floorを利用する．
  floorは必要！もし直接intキャストを使うと，-0.**が0になってしまい，isInsideがfalseなはずがtrueが返ってしまう．

   0.**   1.**   2.**   3.**
  <-dL-> <-dL-> <-dL-> <-dL->
  *-----*------*------*------*
  |  0  |   1  |   2  |   3  |
  *-----*------*------*------*

  */

  // x座標を内包するバケツのインデックスを返す (optimized with cached bounds)
#if defined(__GNUC__) || defined(__clang__)
  [[nodiscard]] __attribute__((always_inline)) inline Tddd itox(const ST i, const ST j, const ST k) const
#else
  [[nodiscard]] inline Tddd itox(const ST i, const ST j, const ST k) const
#endif
  {
    return {this->dL * (0.5 + i) + this->xmin,
            this->dL * (0.5 + j) + this->ymin,
            this->dL * (0.5 + k) + this->zmin};
  }
  [[nodiscard]] inline Tddd itox(const ST3& ijk) const {
    return itox(std::get<0>(ijk), std::get<1>(ijk), std::get<2>(ijk));
  }
  // Optimized: use cached bounds and multiplication instead of division (multiplication is faster)
#if defined(__GNUC__) || defined(__clang__)
  [[nodiscard]] __attribute__((always_inline)) inline ST3 indices_no_clamp(const Tddd& x) const
#else
  [[nodiscard]] inline ST3 indices_no_clamp(const Tddd& x) const
#endif
  {
    //! floorは必要！もし直接intキャストを使うと，-0.**が0になってしまい，isInsideがfalseなはずがtrueが返ってしまう．
    return {static_cast<ST>(std::floor((std::get<0>(x) - this->xmin) * this->inv_dL)),
            static_cast<ST>(std::floor((std::get<1>(x) - this->ymin) * this->inv_dL)),
            static_cast<ST>(std::floor((std::get<2>(x) - this->zmin) * this->inv_dL))};
  };

  [[nodiscard]] inline ST3 indices(const Tddd& x) const {
    const auto ijk = indices_no_clamp(x);
    return {std::clamp(std::get<0>(ijk), static_cast<ST>(0), this->xsize - 1),
            std::clamp(std::get<1>(ijk), static_cast<ST>(0), this->ysize - 1),
            std::clamp(std::get<2>(ijk), static_cast<ST>(0), this->zsize - 1)};
  };
  constexpr ST6 indices_ranges(const Tddd& x, const double d) const {
    auto [i_min, j_min, k_min] = indices(x - d);
    auto [i_max, j_max, k_max] = indices(x + d);
    return {i_min, i_max, j_min, j_max, k_min, k_max};
  };

  constexpr ST6 indices_ranges(const Tddd& x, const Tddd& dx_dy_dz) const {
    auto [dx, dy, dz] = dx_dy_dz;
    auto minmax_x = indices_ranges(x, dx);
    auto minmax_y = indices_ranges(x, dy);
    auto minmax_z = indices_ranges(x, dz);
    return {std::get<0>(minmax_x), std::get<1>(minmax_x), std::get<2>(minmax_y), std::get<3>(minmax_y), std::get<4>(minmax_z), std::get<5>(minmax_z)};
  };

  //% -------------------------------------------------------------------------- */
  //%                                 整合性の確認                                 */
  //% -------------------------------------------------------------------------- */

  // Optimized InsideQ with inline and nodiscard attributes
  // Reordered conditions: check likely-to-fail conditions first (branch prediction)
#if defined(__GNUC__) || defined(__clang__)
  [[nodiscard]] __attribute__((always_inline)) inline bool InsideQ(const ST i, const ST j, const ST k) const
#else
  [[nodiscard]] inline bool InsideQ(const ST i, const ST j, const ST k) const
#endif
  {
    return (i >= 0 && i < this->xsize && j >= 0 && j < this->ysize && k >= 0 && k < this->zsize);
  }

#if defined(__GNUC__) || defined(__clang__)
  [[nodiscard]] __attribute__((always_inline)) inline bool InsideQ(const ST3& ijk) const
#else
  [[nodiscard]] inline bool InsideQ(const ST3& ijk) const
#endif
  {
    return InsideQ(std::get<0>(ijk), std::get<1>(ijk), std::get<2>(ijk));
  }

#if defined(__GNUC__) || defined(__clang__)
  [[nodiscard]] __attribute__((always_inline)) inline bool InsideQ(const Tddd& x) const
#else
  [[nodiscard]] inline bool InsideQ(const Tddd& x) const
#endif
  {
    return InsideQ(indices_no_clamp(x));
  }

  //% -------------------------------------------------------------------------- */

  // Optimized with cached bounds
  [[nodiscard]] inline auto getBounds(const ST3& ijk) const {
    const auto [i, j, k] = ijk;
    return CoordinateBounds(T3Tdd{{{this->dL * i + this->xmin, this->dL * (i + 1) + this->xmin},
                                   {this->dL * j + this->ymin, this->dL * (j + 1) + this->ymin},
                                   {this->dL * k + this->zmin, this->dL * (k + 1) + this->zmin}}});
  };

  //@ -------------------------------------------------------------------------- */

  // unordered_setをvectorにコピー
  void setVector() {
    this->data_vector.assign(this->xsize, std::vector<std::vector<std::vector<T>>>(this->ysize, std::vector<std::vector<T>>(this->zsize)));

    // Iterate through the 3D array and populate the 3D vector
    ST i, j, k;
    for (i = 0; i < this->xsize; ++i) {
      for (j = 0; j < this->ysize; ++j) {
        for (k = 0; k < this->zsize; ++k) {
          auto& bucket = this->data[i][j][k];
          auto& target_vector = this->data_vector[i][j][k];
          target_vector.reserve(bucket.size());
          target_vector.assign(bucket.begin(), bucket.end());
        }
      }
    }

    std::cout << "      bounds : " << this->bounds << std::endl;
    std::cout << "          dL : " << this->dL << std::endl;
    std::cout << "bucket sizes : {" << this->xsize << ", " << this->ysize << ", " << this->zsize << "}" << std::endl;
    this->vector_is_set = true;
  };

  // -------------------------------------------------------------------------- */
  //    erase, shrink, grow, add
  // -------------------------------------------------------------------------- */

  //% erase_add_grow_shrinkは，
  //% deepest_level_bucketsをトラバースしている際に，deepest_level_bucketsが変更になる危険があるため，
  //% それを避けたい場合，erase_addだけを行い，最後にshrinkとgrowを行う．

  //% ------------------------------------------------ */

  void clearDataKeepTree() {
    this->vector_is_set = false;
    this->data1D.clear();
    this->data1D_vector.clear();
    this->map_to_ijk.clear();
    // ツリー構造を保持するためにループしてクリア
    for (auto& plane : this->data)
      for (auto& row : plane)
        for (auto& cell : row)
          cell.clear();
    for (auto& plane : this->data_vector)
      for (auto& row : plane)
        for (auto& cell : row)
          cell.clear();
    traverseTree([](auto& child) { child->clearDataKeepTree(); });
  }

  //% ----------------------------------------------- */

  bool erase_add(const Tddd& x, const T& p) {
    bool success = erase(p);
    success = success && add(x, p);
    return success;
  }

  //% ----------------------------------------------- */

  void shrink() {
    traverseTree([](auto& child) {
      if (child->data1D.empty()) {
        delete child;
        child = nullptr;
      }
    });
  }

  //% ----------------------------------------------- */

  //! 既存subtreeの構造的整合を回復する。
  //! hasChildren()==true の bucket で、非空セルに child がない場合に child を補完する。
  //! grow_condition とは無関係に、既にsplitされたsubtreeを再帰的に閉じる。
  bool repairPartialSubtree() {
    if (!this->hasChildren())
      return false;

    bool changed = false;

    for (int i = 0; i < this->xsize; ++i)
      for (int j = 0; j < this->ysize; ++j)
        for (int k = 0; k < this->zsize; ++k) {
          if (this->data[i][j][k].empty()) {
            // 空セルの子は不要 — shrink で消すべきだが念のため
            if (this->children[i][j][k] != nullptr) {
              delete this->children[i][j][k];
              this->children[i][j][k] = nullptr;
              changed = true;
            }
          } else if (this->children[i][j][k] == nullptr) {
            // 非空セルなのに child がない → 補完
            auto* child = new Buckets<T, N>(getBounds({i, j, k}), this->dL * (0.5 + 1e-13));
            child->setGrowCondition(this->grow_condition);
            child->setLevel(this->level + 1, this->max_level);
            child->parent = this;
            for (auto& p : this->data[i][j][k]) {
              if (!child->add(p->X, p))
                child->add_bypass_insideQ(p->X, p);
            }
            this->children[i][j][k] = child;
            changed = true;
          }
        }

    // 再帰: 子が既にsubtreeを持っている場合、またはgrow_conditionを満たす場合に修復を続行
    for (auto& vi : this->children)
      for (auto& vj : vi)
        for (auto& vk : vj)
          if (vk != nullptr) {
            if (vk->hasChildren() || (vk->grow_condition(vk) && vk->level < vk->max_level))
              changed = vk->repairPartialSubtree() || changed;
          }

    return changed;
  }

  //% ----------------------------------------------- */

  bool grow(bool recursive = true) {
    bool changed = false;

    // Phase 1: 既存subtreeの穴埋め（grow_conditionとは無関係）
    if (this->hasChildren())
      changed = this->repairPartialSubtree();

    // Phase 2: 新規split（grow_conditionに基づく）
    if (!this->grow_condition(this) || this->level >= this->max_level) {
      if (changed && this->level == 0)
        this->rebuildHierarchyLists();
      return changed;
    }

    // 子配列確保 (一度だけ)
    if (!this->hasChildren())
      children.assign(this->xsize, std::vector<std::vector<Buckets<T, N>*>>(this->ysize, std::vector<Buckets<T, N>*>(this->zsize, nullptr)));

    std::vector<std::array<int, 3>> ijk_list;
    ijk_list.reserve(this->xsize * this->ysize * this->zsize);
    for (int i = 0; i < this->xsize; ++i)
      for (int j = 0; j < this->ysize; ++j)
        for (int k = 0; k < this->zsize; ++k)
          ijk_list.emplace_back(std::array<int, 3>{i, j, k});

    // 未生成セルに対して child を生成
    bool local_changed = false;
#pragma omp parallel for reduction(|| : local_changed)
    for (const auto& [i, j, k] : ijk_list) {
      if (this->data[i][j][k].empty()) {
        if (this->children[i][j][k] != nullptr) {
          delete this->children[i][j][k]; // 既存の子は削除
          this->children[i][j][k] = nullptr;
          local_changed = true;
        }
      } else if (children[i][j][k] == nullptr) {
        auto* child = new Buckets<T, N>(getBounds({i, j, k}), this->dL * (0.5 + 1e-13));
        child->setGrowCondition(this->grow_condition);
        child->setLevel(this->level + 1, this->max_level);
        child->parent = this;
        // データ移送 (再帰成長は後段)
        for (auto& p : this->data[i][j][k]) {
          if (!child->add(p->X, p))
            child->add_bypass_insideQ(p->X, p);
        }
        children[i][j][k] = child;
        local_changed = true;
      }
    }
    changed = changed || local_changed;

    // 再帰: 直接の子のみ走査（各子の grow(true) が更に深く再帰する）
    // traverseTree だと全 descendant を回し、各 child->grow(true) も再帰するため重複訪問になる
    if (recursive) {
      traverseChildren([&](Buckets<T, N>*& child) {
        if (child != nullptr && (child->hasChildren() || (child->grow_condition(child) && child->level < child->max_level)))
          changed = child->grow(true) || changed;
      });
    }

    // ルートなら階層リストを再構築
    if (changed && this->level == 0)
      this->rebuildHierarchyLists();

    return changed;
  }

  //% ----------------------------------------------- */

  //% 追加 + 縮小　shrinkを最後に行うことで若干効率が良くなる場合がある．
  //% つまり，入れ直し
  bool erase_add_grow_shrink(const Tddd& x, const T& p) {
    bool success = erase(p);
    success = success && add_grow(x, p);
    traverseTree([](auto& child) {
      if (child->data1D.empty()) {
        delete child;
        child = nullptr;
      }
    });
    //% 入れ直し成功
    return success; //% eraseもaddも成功した場合のみtrueを返す．
  }

  //! 再帰的削除
  bool erase(const T p) {
    bool found = false;
    auto it = this->map_to_ijk.find(p);
    if (it != this->map_to_ijk.end()) {
      this->vector_is_set = false;
      auto [i, j, k] = it->second;
      this->data1D.erase(p);
      auto p_it = std::find(this->data1D_vector.begin(), this->data1D_vector.end(), p);
      if (p_it != this->data1D_vector.end())
        this->data1D_vector.erase(p_it);
      this->data[i][j][k].erase(p);
      this->map_to_ijk.erase(it);
      found = true;
    }
    //! 再帰的削除
    traverseTree([p](auto child) { child->erase(p); });
    return found;
  };

  //! 再帰的削除 + ツリーの縮小
  bool erase_shrink(const T p) {
    const bool found = erase(p);
    traverseTree([](auto& child) {
      if (child->data1D.empty()) {
        delete child;
        child = nullptr;
      }
    });
    return found;
  }

  // 再帰的追加
  bool add(const Tddd& x, const T p) {
    if (!InsideQ(x))
      return false;
    this->vector_is_set = false;
    const auto ijk = this->indices(x);
    this->map_to_ijk[p] = ijk;
    const auto [i, j, k] = ijk;
    const bool bucket_inserted = this->data[i][j][k].emplace(p).second;
    const bool all_objects_inserted = this->data1D.emplace(p).second;
    if (all_objects_inserted)
      this->data1D_vector.emplace_back(p);

    // 親のindicesで決定されたセルの子に直接追加（InsideQ浮動小数点誤差を回避）
    if (this->hasChildren() && children[i][j][k] != nullptr) {
      if (!children[i][j][k]->add(x, p))
        children[i][j][k]->add_bypass_insideQ(x, p);
    }

    return bucket_inserted && all_objects_inserted;
  };

  bool add(const T p) { return this->add(p->X, p); };

  // InsideQチェックをスキップして追加（浮動小数点で境界判定が失敗したソース用）
  void add_bypass_insideQ(const Tddd& x, const T p) {
    this->vector_is_set = false;
    const auto ijk = this->indices(x); // clamp付き → 必ず有効なセルに
    this->map_to_ijk[p] = ijk;
    const auto [i, j, k] = ijk;
    this->data[i][j][k].emplace(p);
    if (this->data1D.emplace(p).second)
      this->data1D_vector.emplace_back(p);
    // 子があれば再帰
    if (this->hasChildren() && children[i][j][k] != nullptr) {
      if (!children[i][j][k]->add(x, p))
        children[i][j][k]->add_bypass_insideQ(x, p);
    }
  };

  // 再帰的追加 + 再帰的追加ツリー生成
  bool add_grow(const Tddd& x, const T& p) {
    if (!InsideQ(x))
      return false;
    const bool ret = this->add(x, p);
    // ツリーがないが，成長条件を満たす場合はツリーを生成し，addする．addしたのでdata1Dがゼロであることはない．
    if (!this->hasChildren() && this->grow_condition(this) && this->level < this->max_level) {
      children.assign(this->xsize, std::vector<std::vector<Buckets<T, N>*>>(this->ysize, std::vector<Buckets<T, N>*>(this->zsize, nullptr)));
      for (int i = 0; i < this->xsize; ++i)
        for (int j = 0; j < this->ysize; ++j)
          for (int k = 0; k < this->zsize; ++k) {
            if (this->data[i][j][k].empty())
              continue;
            Buckets<T, N>* bucket = new Buckets<T, N>(getBounds({i, j, k}), this->dL * (0.5 + 1e-13));
            bucket->setGrowCondition(this->grow_condition); //! 1. Important
            bucket->setLevel(this->level + 1, this->max_level);
            bucket->parent = this; // 親バケツを設定
            children[i][j][k] = bucket;
            for (auto& p : this->data[i][j][k])
              bucket->add_grow(p->X, p); //! 2. Important
          }
    }

    return ret;
  }

  std::tuple<bool, std::array<int, 3>> addAndGetIndices(const Tddd& x, const T p) {
    auto ijk = indices(x);
    return {add(x, p), ijk};
  };

  bool add(const std::unordered_set<T>& P) {
    this->vector_is_set = false;
    bool ret = true;
    for (const auto p : P) {
      if (!add(ToX(p), p))
        ret = false;
    }
    return ret;
  };

  bool add(const std::unordered_set<T>& P, const std::function<Tddd(const T&)>& ToX) {
    this->vector_is_set = false;
    bool ret = true;
    for (const auto p : P) {
      if (!add(ToX(p), p))
        ret = false;
    }
    return ret;
  };

  void initialize_add(const auto& boundingboxIN, const double dL_IN, const std::unordered_set<T>& P) {

    this->initialize(boundingboxIN, dL_IN);
    this->vector_is_set = false;

    for (const auto& p : P) {
      const auto x = ToX(p);
      if (InsideQ(x)) {
        const auto ijk = this->indices(x);
        this->map_to_ijk.emplace(p, ijk); // 事前条件: 未登録なので emplace で十分
        this->data_vector[std::get<0>(ijk)][std::get<1>(ijk)][std::get<2>(ijk)].emplace_back(p);
        this->data[std::get<0>(ijk)][std::get<1>(ijk)][std::get<2>(ijk)].emplace(p);
        traverseTree([x, p](auto child) {
          if (child->InsideQ(x))
            child->add(x, p);
        });
      }
    }

    this->data1D = P;
    this->data1D_vector = std::vector<T>(P.begin(), P.end());
    this->vector_is_set = true;
  };

  /* -------------------------------------------------------------------------- */

  auto getBucket(const Tddd& x) const {
    const auto ijk = this->indices(x);
    return this->children[std::get<0>(ijk)][std::get<1>(ijk)][std::get<2>(ijk)];
  };

  auto getData(const Tddd& x) const {
    const auto ijk = this->indices(x);
    return this->data[std::get<0>(ijk)][std::get<1>(ijk)][std::get<2>(ijk)];
  };

  auto getData(const Tddd& x, const double& range) const {
    const auto ijk = this->indices(x);
    std::unordered_set<T> ret;
    this->apply(x, range, [&](auto& a) { ret.emplace(a); });
    return ret;
  };

  void clear() {
    this->vector_is_set = false;
    this->data.clear();
    this->data_vector.clear();
    this->data1D.clear();
    this->data1D_vector.clear();
    this->map_to_ijk.clear();
    this->is_data_assigned_as_3D = false; // 初期状態に戻すこと
  };

  // b@ ========================================================================= */
  // b@　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　localなdataに対する　STL like　functions　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　*/
  // b@ ========================================================================= */

  //! none_of
  bool none_of(const Tddd& x, const double d, const std::function<bool(const T&)>& func) const {
    if (this->data.empty())
      return true;
    const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);
    if (!this->vector_is_set) {
      return !spatial_partitioning_detail::any_of_unseq(this->data.cbegin() + i_min, this->data.cbegin() + i_max + 1, [&](const auto& Bi) { return spatial_partitioning_detail::any_of_unseq(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto& Bij) { return spatial_partitioning_detail::any_of_unseq(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto& Bijk) { return spatial_partitioning_detail::any_of_unseq(Bijk.cbegin(), Bijk.cend(), [&](const auto p) { return func(p); }); }); }); });
    } else {
      return !spatial_partitioning_detail::any_of_unseq(this->data_vector.cbegin() + i_min, this->data_vector.cbegin() + i_max + 1, [&](const auto& Bi) { return spatial_partitioning_detail::any_of_unseq(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto& Bij) { return spatial_partitioning_detail::any_of_unseq(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto& Bijk) { return spatial_partitioning_detail::any_of_unseq(Bijk.cbegin(), Bijk.cend(), [&](const auto p) { return func(p); }); }); }); });
    }
  }

  //! all_of
  bool all_of(const Tddd& x, const double d, const std::function<bool(const T&)>& func) const {
    if (this->data.empty())
      return true;
    const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);
    if (!this->vector_is_set) {
      return spatial_partitioning_detail::all_of_unseq(this->data.cbegin() + i_min, this->data.cbegin() + i_max + 1, [&](const auto& Bi) { return spatial_partitioning_detail::all_of_unseq(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto& Bij) { return spatial_partitioning_detail::all_of_unseq(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto& Bijk) { return spatial_partitioning_detail::all_of_unseq(Bijk.cbegin(), Bijk.cend(), [&](const auto p) { return func(p); }); }); }); });
    } else {
      return spatial_partitioning_detail::all_of_unseq(this->data_vector.cbegin() + i_min, this->data_vector.cbegin() + i_max + 1, [&](const auto& Bi) { return spatial_partitioning_detail::all_of_unseq(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto& Bij) { return spatial_partitioning_detail::all_of_unseq(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto& Bijk) { return spatial_partitioning_detail::all_of_unseq(Bijk.cbegin(), Bijk.cend(), [&](const auto p) { return func(p); }); }); }); });
    }
  }

  //! any_of
  bool any_of(const Tddd& x, const double d, const std::function<bool(const T&)>& func) const {
    if (this->data.empty())
      return false;
    const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);
    if (!this->vector_is_set) {
      return spatial_partitioning_detail::any_of_unseq(this->data.cbegin() + i_min, this->data.cbegin() + i_max + 1, [&](const auto& Bi) { return spatial_partitioning_detail::any_of_unseq(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto& Bij) { return spatial_partitioning_detail::any_of_unseq(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto& Bijk) { return spatial_partitioning_detail::any_of_unseq(Bijk.cbegin(), Bijk.cend(), [&](const auto p) { return func(p); }); }); }); });
    } else {
      return spatial_partitioning_detail::any_of_unseq(this->data_vector.cbegin() + i_min, this->data_vector.cbegin() + i_max + 1, [&](const auto& Bi) { return spatial_partitioning_detail::any_of_unseq(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto& Bij) { return spatial_partitioning_detail::any_of_unseq(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto& Bijk) { return spatial_partitioning_detail::any_of_unseq(Bijk.cbegin(), Bijk.cend(), [&](const auto p) { return func(p); }); }); }); });
    }
  }

  //! apply
  void apply(const Tddd& x, const double d, const std::function<void(const T&)>& func) const {
    if (this->data.empty())
      return;
    const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);
    if (!this->vector_is_set) {
      spatial_partitioning_detail::for_each_unseq(this->data.cbegin() + i_min, this->data.cbegin() + i_max + 1,
                                                  [&func, &j_min, &j_max, &k_min, &k_max](const auto& Bi) { spatial_partitioning_detail::for_each_unseq(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&func, &k_min, &k_max](const auto& Bij) { spatial_partitioning_detail::for_each_unseq(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&func](const auto& Bijk) { spatial_partitioning_detail::for_each_unseq(Bijk.cbegin(), Bijk.cend(), [&](const auto p) { func(p); }); }); }); });
    } else {
      spatial_partitioning_detail::for_each_unseq(this->data_vector.cbegin() + i_min, this->data_vector.cbegin() + i_max + 1,
                                                  [&func, &j_min, &j_max, &k_min, &k_max](const auto& Bi) { spatial_partitioning_detail::for_each_unseq(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&func, &k_min, &k_max](const auto& Bij) { spatial_partitioning_detail::for_each_unseq(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&func](const auto& Bijk) { spatial_partitioning_detail::for_each_unseq(Bijk.cbegin(), Bijk.cend(), [&](const auto p) { func(p); }); }); }); });
    }
  }

  //! apply
  void apply(const Tddd& x, const double d, const std::function<void(const int, const int, const int)>& func) const {
    const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);
    int i, j, k;
    for (i = i_min; i <= i_max; ++i)
      for (j = j_min; j <= j_max; ++j)
        for (k = k_min; k <= k_max; ++k)
          func(i, j, k);
  }

  void applyVec1D(const Tddd& x, const double d, const std::function<void(const std::vector<T>&)>& func) const {
    if (this->data.empty() || !this->vector_is_set)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "'s 3D data is empty");
    const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);
    spatial_partitioning_detail::for_each_unseq(this->data_vector.cbegin() + i_min, this->data_vector.cbegin() + i_max + 1, [&func, &j_min, &j_max, &k_min, &k_max](const auto& Bi) { spatial_partitioning_detail::for_each_unseq(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&func, &k_min, &k_max](const auto& Bij) { spatial_partitioning_detail::for_each_unseq(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&func](const auto& Bijk) { func(Bijk); }); }); });
  }

  //! apply
  void apply_to_the_nearest_bound(const Tddd& x, const std::function<void(const int, const int, const int)>& func) const {
    if (this->data.empty())
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "'s 3D data is empty");
    auto [I, J, K] = indices(x);

    int values[] = {I, J, K, this->xsize - I, this->ysize - J, this->zsize - K};
    int smallest = *std::min_element(values, values + 6);
    int smallestIndex = std::distance(values, std::find(values, values + 6, smallest));

    int i, j, k;

    if (smallestIndex == 0) {
      for (i = I; i >= 0; --i)
        func(i, J, K);
    } else if (smallestIndex == 1) {
      for (j = J; j >= 0; --j)
        func(I, j, K);
    } else if (smallestIndex == 2) {
      for (k = K; k >= 0; --k)
        func(I, J, k);
    } else if (smallestIndex == 3) {
      for (i = I; i < this->xsize; ++i)
        func(i, J, K);
    } else if (smallestIndex == 4) {
      for (j = J; j < this->ysize; ++j)
        func(I, j, K);
    } else if (smallestIndex == 5) {
      for (k = K; k < this->zsize; ++k)
        func(I, J, k);
    }
  }

  //! apply
  void apply(const Tddd& x, const Tddd dx_dy_dz, const std::function<void(const T&)>& func) const {
    if (this->data.empty())
      return;
    const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, dx_dy_dz);
    if (!this->vector_is_set) {
      spatial_partitioning_detail::for_each_unseq(this->data.cbegin() + i_min, this->data.cbegin() + i_max + 1,
                                                  [&func, &j_min, &j_max, &k_min, &k_max](const auto& Bi) { spatial_partitioning_detail::for_each_unseq(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&func, &k_min, &k_max](const auto& Bij) { spatial_partitioning_detail::for_each_unseq(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&func](const auto& Bijk) { spatial_partitioning_detail::for_each_unseq(Bijk.cbegin(), Bijk.cend(), [&](const auto p) { func(p); }); }); }); });
    } else {
      spatial_partitioning_detail::for_each_unseq(this->data_vector.cbegin() + i_min, this->data_vector.cbegin() + i_max + 1,
                                                  [&func, &j_min, &j_max, &k_min, &k_max](const auto& Bi) { spatial_partitioning_detail::for_each_unseq(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&func, &k_min, &k_max](const auto& Bij) { spatial_partitioning_detail::for_each_unseq(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&func](const auto& Bijk) { spatial_partitioning_detail::for_each_unseq(Bijk.cbegin(), Bijk.cend(), [&](const auto p) { func(p); }); }); }); });
    }
  };

  /* ----------------------------- FOR LINE SEARCH ---------------------------- */

  std::vector<ST3> line2indices(const Tddd& A, const Tddd& B) const {
    std::unordered_set<ST3> uniqueIndexSet;
    for (const auto& X : Subdivide(A, B, std::ceil(Norm(A - B) / (this->dL / 2.))))
      this->apply(X, this->dL, [&uniqueIndexSet](const int i, const int j, const int k) { uniqueIndexSet.insert({i, j, k}); });
    std::vector<ST3> uniqueIndices(uniqueIndexSet.begin(), uniqueIndexSet.end());
    return uniqueIndices;
  }

  //! apply using indices
  void apply(const std::vector<ST3>& V_ijk, const std::function<void(const T&)>& func) const {
    if (!this->vector_is_set) {
      for (const auto& ijk : V_ijk)
        for (const auto& p : this->data[std::get<0>(ijk)][std::get<1>(ijk)][std::get<2>(ijk)])
          func(p);
    } else {
      for (const auto& ijk : V_ijk)
        for (const auto& p : this->data_vector[std::get<0>(ijk)][std::get<1>(ijk)][std::get<2>(ijk)])
          func(p);
    }
  }

  //^ ============================================================================　*/
  //^　　　　　　　　　　　　　　　　　　　ツリー構造　　　　　　　　　　　　　　　　　　　　　　　*/
  //^　　　　　　　　childrenはすでにあるdataの分割を真似るようにして作成される　　　　　　　　　*/
  //^ ============================================================================　*/

  std::vector<std::vector<std::vector<Buckets<T, N>*>>> children;
  Buckets<T, N>* parent = nullptr; // 親バケツへのポインタを追加

  std::vector<Buckets<T, N>*> getAllChildren() {
    std::vector<Buckets<T, N>*> all_buckets;
    traverseChildren([&](Buckets<T, N>* b) { all_buckets.emplace_back(b); });
    return all_buckets;
  }

  /* -------------------------------------------------------------------------- */

  bool hasChildren() const {
    for (const auto& vi : this->children)
      for (const auto& vj : vi)
        for (const auto& vk : vj)
          if (vk != nullptr)
            return true;
    return false;
  }
  int level = 0;     //! ツリー全体における，このバケツのレベル
  int max_level = 1; //! ツリーの最深レベル

  void setLevel(const int levelIN, const int max_levelIN) {
    this->level = levelIN;
    this->max_level = max_levelIN;
  }

  /* ---------------------------------------------------- ---------------------- */
  /* -------------------------------------------------------------------------- */
  /* -------------------------------------------------------------------------- */

  std::vector<std::vector<Buckets<T, N>*>> level_buckets;
  std::vector<Buckets<T, N>*> deepest_level_buckets;

  //! 無駄がないようにgenerate Treeしよう！
  void rebuildHierarchyLists() {
    this->level_buckets.clear();
    this->deepest_level_buckets.clear();
    this->level_buckets.resize(this->max_level + 1);
    //
    if (!this->data1D.empty()) {
      this->level_buckets[this->level].emplace_back(this);
      if (!this->hasChildren())
        this->deepest_level_buckets.emplace_back(this);
    }
    //
    this->traverseTree([&](Buckets<T, N>* b) -> void {
      if (!b->data1D.empty()) {
        this->level_buckets[b->level].emplace_back(b);
        if (!b->hasChildren())
          this->deepest_level_buckets.emplace_back(b);
      }
    });
  };

  /* -------------------------------------------------------------------------- */
  /* -------------------------------------------------------------------------- */
  /* -------------------------------------------------------------------------- */

  std::vector<T> rebin() {
    this->vector_is_set = false;
    int success = 0;
    auto tryMigrateUpward = [&success](Buckets<T, N>* target, const std::vector<T>& objectsOutOfBounds, auto& tryMigrate_) -> std::vector<T> {
      if (target == nullptr || objectsOutOfBounds.empty())
        return {};
      //
      std::vector<T> pass2parent;
      pass2parent.reserve(objectsOutOfBounds.size());
      for (auto& p : objectsOutOfBounds) {
        if (!target->erase_add(p->X, p))
          pass2parent.emplace_back(p);
        else
          success++;
      }
      //
      if (target->parent == nullptr || pass2parent.empty())
        return pass2parent; // 親がいない場合はここで終了
      return tryMigrate_(target->parent, pass2parent, tryMigrate_);
    };

    std::vector<T> escaped, escaped_once;
    // Reserve to avoid repeated allocations
    escaped.reserve(this->deepest_level_buckets.size());
    escaped_once.reserve(this->deepest_level_buckets.size());

    for (auto* b : this->deepest_level_buckets) {
      std::vector<T> escapedAtLeaf;
      escapedAtLeaf.reserve(b->data1D.size()); // Reserve based on bucket size
      for (auto p : b->data1D)
        if (!b->InsideQ(p->X)) {
          escapedAtLeaf.emplace_back(p);
          escaped_once.emplace_back(p);
        }
      for (auto& p : tryMigrateUpward(b->parent, escapedAtLeaf, tryMigrateUpward))
        escaped.emplace_back(p);
    }

    std::cout << Green << "Escaped objects after rebin: " << escaped_once.size() << colorReset << std::endl;
    std::cout << Yellow << "Success objects after rebin: " << success << colorReset << std::endl;
    std::cout << Red << "Remaining objects after rebin: " << escaped.size() << colorReset << std::endl;

    // If nothing escaped, there's no structural change; skip expensive shrink/grow.
    if (escaped_once.empty())
      return {};

    // deepest_level_bucketsの更新
    if (this->level == 0) {
      this->rebuildHierarchyLists();
      std::cout << "Before shrink: " << this->deepest_level_buckets.size() << " deepest level buckets." << std::endl;
      this->shrink();
      // IMPORTANT:
      // `shrink()` deletes bucket nodes, so any cached pointers (deepest_level_buckets/level_buckets)
      // may contain dangling pointers unless we rebuild them right away.
      this->rebuildHierarchyLists();
      bool changed = this->grow();
      // `grow()` rebuilds hierarchy lists only when it creates new nodes.
      // If shrink deleted nodes but grow doesn't create any, we still need a rebuild.
      if (changed)
        this->rebuildHierarchyLists();
      std::cout << "After grow: " << this->deepest_level_buckets.size() << " deepest level buckets." << std::endl;
    }
    // rebin後のリーフソース数検証 — mismatch時はfull rebuildにfallback
    {
      std::size_t total_leaf_sources = 0;
      for (auto* b : this->deepest_level_buckets)
        total_leaf_sources += b->data1D.size();
      if (total_leaf_sources != this->data1D.size()) {
        std::cerr << Red << "[rebin] WARNING: leaf sources (" << total_leaf_sources
                  << ") != total sources (" << this->data1D.size()
                  << "). Falling back to full rebuild." << colorReset << std::endl;
        // 子ツリーを全削除（rootのdata/data1Dは保持）
        for (auto& vi : this->children)
          for (auto& vj : vi)
            for (auto& vk : vj)
              if (vk != nullptr) {
                delete vk;
                vk = nullptr;
              }
        this->children.clear();
        // rootのdata1Dから木を再構築
        this->generateTree();
      }
    }
    return escaped;
  }

  /* -------------------------------------------------------------------------- */
  /* -------------------------------------------------------------------------- */
  /* -------------------------------------------------------------------------- */

private:
  std::function<bool(const Buckets<T, N>*)> default_grow_condition = [](const Buckets<T, N>* b) { return !b->data1D.empty() && b->level < b->max_level; };
  std::function<bool(const Buckets<T, N>*)> grow_condition = default_grow_condition;
  bool grow_condition_setQ = false;

public:
  void setGrowCondition(const std::function<bool(const Buckets<T, N>*)>& condition) {
    this->grow_condition = condition;
    this->grow_condition_setQ = true;
  }

  bool generateTree() {
    if (this->grow_condition_setQ)
      return this->generateTree(this->grow_condition);
    else
      return this->generateTree(this->default_grow_condition);
  }

private:
  Buckets<T, N>* newChild(int i, int j, int k) {
    Buckets<T, N>* bucket = new Buckets<T, N>(getBounds({i, j, k}), this->dL * (0.5 + 1e-13));
    bucket->parent = this; // 親バケツを設定
    for (auto& p : this->data[i][j][k])
      bucket->add(p->X, p);
    bucket->setLevel(this->level + 1, this->max_level);
    return bucket;
  }

  bool generateTree(const std::function<bool(const Buckets<T, N>*)>& condition) {
    this->grow_condition = condition; // Automatically update grow_condition
    if (condition(this) && this->level < this->max_level) {
      children.assign(this->xsize, std::vector<std::vector<Buckets<T, N>*>>(this->ysize, std::vector<Buckets<T, N>*>(this->zsize, nullptr)));
#pragma omp parallel for schedule(dynamic) if (this->level < 5)
      for (int i = 0; i < (int)this->vector_ijk.size(); ++i) {
        auto [ix, iy, iz] = this->vector_ijk[i];
        if (this->data[ix][iy][iz].empty())
          continue;
        children[ix][iy][iz] = newChild(ix, iy, iz);
        children[ix][iy][iz]->generateTree(condition);
      }
    }

    if (this->level == 0)
      rebuildHierarchyLists();
    return this->hasChildren();
  }

  //! ツリー構造を維持したまま、データ点の移動に合わせてツリーを更新する
  //! generateTreeの代わりにこれを呼ぶことで高速化が期待できる
  void updateTree() {
    if (this->level != 0) {
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "updateTree must be called from the root bucket.");
    }
    // 移動した粒子を適切なバケットに移動させ、空になった枝を削除し、必要に応じて枝を伸ばす
    [[maybe_unused]] auto escaped = this->rebin();
  }

public:
  /* =================== Access to the bucket at the level of the tree =================== */

  void traverseChildren(const std::function<void(Buckets<T, N>*&)>& func) {
    for (auto& i : this->children)
      for (auto& j : i)
        for (auto& k : j)
          if (k != nullptr)
            func(k);
  }

  void traverseTree(const std::function<void(Buckets<T, N>*&)>& func) {
    for (auto& i : this->children)
      for (auto& j : i)
        for (auto& k : j)
          if (k != nullptr) {
            func(k);
            if (k != nullptr)
              k->traverseTree(func);
          }
  }

  void traverseTree(const std::function<void(Buckets<T, N>*&)>& func, const std::function<bool(Buckets<T, N>*&)>& stop_condition) {
    for (auto& i : this->children)
      for (auto& j : i)
        for (auto& k : j)
          if (k != nullptr && !stop_condition(k)) {
            func(k);
            if (k != nullptr)
              k->traverseTree(func, stop_condition);
          }
  }

  void traverseTreeParallel(const std::function<void(Buckets<T, N>*&)>& func) {
    spatial_partitioning_detail::for_each_par_unseq(this->children.begin(), this->children.end(), [&](auto& i) {
      spatial_partitioning_detail::for_each_unseq(i.begin(), i.end(), [&](auto& j) {
        spatial_partitioning_detail::for_each_unseq(j.begin(), j.end(), [&](auto& k) {
          if (k != nullptr) {
            func(k);
            if (k != nullptr)
              k->traverseTree(func);
          }
        });
      });
    });
  }

  //@ applyToBucket(1,func)
  //@ bucketsは，level==1のバケツに適用される

  //@ this is traversing using the tree structure.これにlevel_bucketsを使ってはならない．なぜなら，level_bucketsの生成にこの関数が使われるから．
  void forEachAll(const std::function<void(Buckets<T, N>*)>& func_for_bucket) {
    func_for_bucket(this);
    spatial_partitioning_detail::for_each_unseq(this->children.begin(), this->children.end(), [&](auto& i) {
      spatial_partitioning_detail::for_each_unseq(i.begin(), i.end(), [&](auto& j) {
        spatial_partitioning_detail::for_each_unseq(j.begin(), j.end(), [&](auto& k) {
          if (k != nullptr) {
            k->forEachAll(func_for_bucket);
          }
        });
      });
    });
  }

  void forEachAllParallel(const std::function<void(Buckets<T, N>*)>& func_for_bucket) {
    func_for_bucket(this);
    spatial_partitioning_detail::for_each_par_unseq(this->children.begin(), this->children.end(), [&](auto& i) {
      spatial_partitioning_detail::for_each_unseq(i.begin(), i.end(), [&](auto& j) {
        spatial_partitioning_detail::for_each_unseq(j.begin(), j.end(), [&](auto& k) {
          if (k != nullptr) {
            k->forEachAllParallel(func_for_bucket);
          }
        });
      });
    });
  }

  void forEachAtLevel(const int level, const std::function<void(Buckets<T, N>*)>& func_for_bucket) {
    if (level >= this->level_buckets.size())
      return;

    for (auto& b : this->level_buckets[level])
      func_for_bucket(b);
  }

  void forEachAtLevel(const std::vector<int>& levels, const std::function<void(Buckets<T, N>*)>& func_for_bucket) {
    for (const auto& lv : levels)
      forEachAtLevel(lv, func_for_bucket);
  }

  void forEachAtLevelParallel(const int lv, const std::function<void(Buckets<T, N>*)>& func_for_bucket) {
    if (lv >= this->level_buckets.size())
      return;
#pragma omp parallel for
    for (auto& b : this->level_buckets[lv])
      func_for_bucket(b);
  }

  void forEachAtLevelParallel(const std::vector<int>& levels, const std::function<void(Buckets<T, N>*)>& func_for_bucket) {
    //! this->level_bucketsは，{0,1,2,3}のように，levelに対応するバケツのリストを持っている．(size()=4)
    //! この場合，level=4を指定するとエラー，これをチェックするには，lv < this->level_buckets.size()とする．
    for (const auto& lv : levels)
      forEachAtLevelParallel(lv, func_for_bucket);
  }

  template <typename Func>
  void forEachAtDeepestParallel(const Func& func_for_bucket) {
#pragma omp parallel for
    for (const auto& b : this->deepest_level_buckets) {
      {
        func_for_bucket(b);
      }
    }
  }

  template <typename Func>
  void forEachAtDeepest(const Func& func_for_bucket) {
    for (auto& b : this->deepest_level_buckets)
      func_for_bucket(b);
  }

  void forEachBuckets(const std::function<void(Buckets<T, N>*)>& func_for_bucket) {
    for (auto& i : this->children)
      for (auto& j : i)
        for (auto& k : j)
          if (k != nullptr)
            func_for_bucket(k);
  }

  Buckets<T, N>* getBucketAtLevel(const int level, const Tddd& x) const {
    if (!this->hasChildren() || level > this->max_level)
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "The tree structure does not exist or the level is too large.");

    auto current = const_cast<Buckets<T, N>*>(this);
    while (current->level < level) {
      if (!current->hasChildren())
        return nullptr;
      auto next = current->getBucket(x);
      if (next == nullptr)
        return nullptr;
      current = next;
    }
    return current;
  };

  Buckets<T, N>* getBucketAtDeepest(const Tddd& x) const {
    auto current = const_cast<Buckets<T, N>*>(this);
    while (current->hasChildren()) {
      auto next = current->getBucket(x);
      if (next == nullptr)
        return current;
      current = next;
    }
    return current;
  };

  // b$============================================================================*/
  // b$                         高速多重極展開に関連するメンバ                           */
  // b$============================================================================*/

  using MomentsType = Moments<T, N>;
  static_assert(std::is_trivially_destructible_v<MomentsType> || std::is_destructible_v<MomentsType>, "Moments must be destructible");
  // インターフェース存在チェック (テンプレート引数を伴うものは適切な引数で検証)
  static_assert(requires(MomentsType m, std::vector<Buckets<T, N>*> v) { m.set_m2m(v); }, "Missing Moments::set_m2m(std::vector<Buckets*>)");
  static_assert(requires(MomentsType m) { m.m2m(); }, "Missing Moments::m2m()");
  //    static_assert(requires(MomentsType m, std::vector<Buckets<T, N> *> v) { m.set_m2l(v); }, "Missing Moments::set_m2l(std::vector<Buckets*>)");
  static_assert(requires(MomentsType m) { m.m2l(); }, "Missing Moments::m2l()");
  static_assert(requires(MomentsType m) { m.set_l2l(&m); }, "Missing Moments::set_l2l(Moments*)");
  static_assert(requires(MomentsType m) { m.l2l(); }, "Missing Moments::l2l()");

  MomentsType MomentsMultipoleExpansion;
  MomentsType MomentsLocalExpansion;

  //$ buckets_for_M2Lは，このバケツが局所展開する対象を保存している．
  std::vector<Buckets<T, N>*> buckets_for_M2L;
  //$ buckets_for_L2Mは，このバケツが局所展開を受け取る対象を保存している．
  std::vector<Buckets<T, N>*> buckets_for_L2M;
  std::vector<Buckets<T, N>*> buckets_near; //@ 近傍バケツかつ，（重要）子バケツを持たないバケツ．途中までしか育たないツリーを使ったFMMにおいてこの性質が重要．;

  // MAC拡張用: バケツ内の全要素の max_source_offset の最大値
  double max_internal_offset = 0.0;

  // b$=============================================================================*/
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

template <typename T>
Buckets<T> copyPartition(const auto& children) {
  Buckets<T> ret(children.bounds, children.dL);
  ret.children.resize(children.xsize, std::vector<std::vector<Buckets<T>*>>(children.ysize, std::vector<Buckets<T>*>(children.zsize, nullptr)));
  ret.level = children.level;
  ret.max_level = children.max_level;
  ret.xsize = children.xsize;
  ret.ysize = children.ysize;
  ret.zsize = children.zsize;
  ret.center = children.center;
  ret.dn = children.dn;
  ret.data.resize(children.xsize, std::vector<std::vector<std::unordered_set<T>>>(children.ysize, std::vector<std::unordered_set<T>>(children.zsize, std::unordered_set<T>{})));
  ret.data_vector.resize(children.xsize, std::vector<std::vector<std::vector<T>>>(children.ysize, std::vector<std::vector<T>>(children.zsize)));
  ret.vector_is_set = false;
  return ret;
}
