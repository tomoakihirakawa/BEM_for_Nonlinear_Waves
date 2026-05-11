#pragma once

#include "Network.hpp"
#include <algorithm>
#include <cstdint>

/*
English terminology
===================

This header is the canonical place where naming is fixed for the dynamic
boundary-condition implementation. The same terminology should be used
consistently in:

  - code comments and new function names,
  - technical notes and documentation,
  - paper drafts and presentations.

Geometry vs. nodes
------------------

  - "vertex" means a geometric corner of a surface triangle.
  - "edge" means a geometric edge of a surface triangle.

  - "vertex node" means a BEM node located on a vertex.
    In the current code this is represented by networkPoint.

  - "edge node" means a BEM node located on an edge.
    In the current code this is represented by networkLine when it is used as
    the midpoint node of a true/pseudo quadratic boundary element.

The word "node" alone is intentionally reserved for the unified concept
"vertex node or edge node". When distinction matters, always say
"vertex node" or "edge node".

Node-face unit
--------------

The primary local unit of boundary classification is not a face by itself.
It is a pair:

  - (vertex node, face)
  - (edge node, face)

In prose, these should be called:

  - "vertex-node face state"
  - "edge-node face state"

or collectively:

  - "node-face state"

Use "state" when the object includes more than a numerical quantity alone, for
example:

  - Dirichlet / Neumann type,
  - pressure-detachment flag,
  - detach counter or eligibility,
  - supporting contact information.

Use "quantity" when the object is a numerical quantity attached to a node-face
pair. A quantity may be scalar, vector, or tensor. For example:

  - phi on a node-face pair,
  - phin on a node-face pair,
  - pressure on a node-face pair,
  - velocity on a node-face pair,
  - stress tensor on a node-face pair.

Use "value" only for the concrete value of a quantity. When needed, say
"scalar value", "vector value", or "tensor value". If the context already
makes it clear which quantity is meant, "vector value" or "tensor value"
alone is acceptable.

Therefore:

  - boundary-condition kind => node-face state
  - phi / phin / pressure   => node-face quantity
  - velocity / stress       => node-face quantity
  - concrete numerical data => value of a node-face quantity

Why this matters
----------------

The implementation allows different boundary information on the same geometric
location depending on which face is considered. Therefore, face/line/point
summary flags are derived quantities only. The primary local information lives
on node-face pairs and should be described with the terminology above.

Japanese terminology
====================

このヘッダを，動的境界条件まわりの用語規約の正本とする．今後は以下の日本語を
コード，技術文書，論文草稿，発表資料で統一して使う．

幾何と節点
----------

  - "vertex"            = 頂点
  - "edge"              = 辺
  - "vertex node"       = 頂点節点
  - "edge node"         = 辺節点
  - "node"              = 節点

ここで，

  - 頂点節点は頂点上に置かれた BEM 節点であり，
    現在のコードでは networkPoint に対応する．

  - 辺節点は辺上に置かれた BEM 節点であり，
    現在のコードでは true/pseudo quadratic 要素の中点節点として使われる
    networkLine に対応する．

"節点" 単独という語は，頂点節点と辺節点をまとめた総称としてのみ使う．
区別が必要なときは，必ず "頂点節点" または "辺節点" と書く．

節点-面単位
-----------

境界条件の一次的な局所単位は face 単独ではなく，次の pair である．

  - (頂点節点, face)
  - (辺節点, face)

本文・コメントでは，これらを次のように呼ぶ．

  - "vertex-node face state" = 頂点節点-面ペアの局所境界状態
  - "edge-node face state"   = 辺節点-面ペアの局所境界状態

総称としては，

  - "node-face state"        = 節点-面状態 / 局所境界状態

を用いる．

"state" と "quantity" と "value" の使い分け
--------------------------------------------

"state" は，単なる量 1 個ではなく，境界条件種別や履歴・フラグも含む局所情報
全体を指す．例えば，

  - Dirichlet / Neumann 種別
  - pressure-detachment flag
  - detach counter や eligibility
  - 接触相手の情報

を含むものは "state" と呼ぶ．

"quantity" は，節点-面ペアに付随する量そのものを指す．"quantity" は
スカラー量・ベクトル量・テンソル量をすべて含む．例えば，

  - 節点-面ペア上の phi
  - 節点-面ペア上の phin
  - 節点-面ペア上の pressure
  - 節点-面ペア上の velocity
  - 節点-面ペア上の stress tensor

は "quantity" と呼ぶ．

"value" は，その quantity の具体的な値を指す．必要に応じて

  - scalar value
  - vector value
  - tensor value

のように明示して使う．文脈上どの quantity の値かが自明な場合には，
"vector value" や "tensor value" とだけ書いてもよい．

したがって，

  - 境界条件の種別や剥離状態     => node-face state
  - phi / phin / pressure        => node-face quantity
  - ベクトル量・テンソル量       => node-face quantity
  - 量の具体的な数値            => value

とする．

この区別が必要な理由
--------------------

この実装では，同じ幾何学的位置であっても，どの face 側を見るかによって異なる
境界情報を持ちうる．したがって，face/line/point の summary flag は派生量に
すぎず，一次的な局所情報は節点-面ペア上にある．用語もこの構造を反映して
記述すること．
*/

// Primary boundary-condition kind carried by a node-face state.
// face/line/point Dirichlet/Neumann summaries are derived from this state.
enum class NodeFaceBoundaryType : std::uint8_t {
  Undefined = 0,
  Dirichlet = 1,
  Neumann = 2,
};

template <class Entity>
inline bool isAdjacentEntityFace(const Entity* entity, const networkFace* f) {
  if (!entity || !f)
    return false;
  const auto bfs = entity->getBoundaryFaces();
  return std::find(bfs.begin(), bfs.end(), f) != bfs.end();
}

template <class Entity>
inline bool isEntityFaceDetachedByPressure(const Entity* entity, const networkFace* f) {
  if (!entity || !f)
    return false;
  const auto* d = entity->findContactState(f);
  return d && d->detached_by_pressure;
}

template <class Entity>
inline NodeFaceBoundaryType getNodeFaceBoundaryType(const Entity* entity, const networkFace* f) {
  if (!isAdjacentEntityFace(entity, f))
    return NodeFaceBoundaryType::Undefined;

  const auto* d = entity->findContactState(f);
  if (d && d->isNeumann())
    return NodeFaceBoundaryType::Neumann;

  // Fallback: if the face centroid is inside a rigid body (penetratedBody),
  // treat as Neumann unless pressure-detachment overrides.
  // This restores the old InsideQ → f->Neumann → p->Neumann path.
  // No circular dependency: updateFacePenetratedBodiesFromTraceStates resets
  // f->penetratedBody = nullptr before calling isNeumannBoundaryState.
  if (f && f->penetratedBody && !isEntityFaceDetachedByPressure(entity, f))
    return NodeFaceBoundaryType::Neumann;

  return NodeFaceBoundaryType::Dirichlet;
}

inline bool isAdjacentPointFace(const networkPoint* p, const networkFace* f) {
  return isAdjacentEntityFace(p, f);
}

inline bool isAdjacentLineFace(const networkLine* l, const networkFace* f) {
  return isAdjacentEntityFace(l, f);
}

inline NodeFaceBoundaryType getVertexNodeFaceBoundaryType(const networkPoint* p, const networkFace* f) {
  return getNodeFaceBoundaryType(p, f);
}

inline NodeFaceBoundaryType getEdgeNodeFaceBoundaryType(const networkLine* l, const networkFace* f) {
  return getNodeFaceBoundaryType(l, f);
}

/* ---------------------------------------------------------------------------
   Layer 1: BoundaryState — 接触状態から決まる境界種別（「状態」）

   境界条件判定の3層構造:

     Layer 1: BoundaryState   — 接触状態としてどうなっているか（このファイル）
     Layer 2: BieDofKey       — BIE unknown の代表キーとしてどの組を使うべきか
                                （BEM_setBoundaryTypes.hpp: isNeumannBieDofKey 等）
     Layer 3: ActiveBieDof    — 実際に index >= 0 の DOF が存在するか
                                （Network.hpp: findActiveBieDof 等）

   利用可能タイミング:
     setContactFaces() 完了後（接触状態が確定した後）。
     setNodeFaceIndices() の前後どちらでも使える。

   Layer 1 は接触検出の結果のみに基づき、DOF の index は参照しない。
   「この (entity, face) ペアは Neumann 境界状態か」を問う。

   point (networkPoint) と line (networkLine) の両方で使えるテンプレート。
   旧名: isVertexNodeFaceNeumann / isEdgeNodeFaceNeumann
         isVertexNodeFaceDirichlet / isEdgeNodeFaceDirichlet
--------------------------------------------------------------------------- */
inline bool isNeumannBoundaryState(const auto* entity, const networkFace* f) {
  return getNodeFaceBoundaryType(entity, f) == NodeFaceBoundaryType::Neumann;
}

inline bool isDirichletBoundaryState(const auto* entity, const networkFace* f) {
  return getNodeFaceBoundaryType(entity, f) == NodeFaceBoundaryType::Dirichlet;
}

/* ---------------------------------------------------------------------------
   Multiple node detection and assignment

   isMultipleNode = true の節点には per-face DOF が割り当てられる。
   判定基準は**SharpQ と同じ幾何 crease**:

     境界辺の隣接面法線同士の最大角が閾値角度を超える → multiple

   BCInterface (D-N接合) であるかどうかは判定に使わない。
   BCInterface は必然的に鋭角（壁面と自由表面の接合）なので、
   法線角度の判定だけで自然に multiple になる。

   ただし、Dirichlet のみの節点は現在の実装では multiple にできない。
   Neumann の multiple node では BIE の式を面ごとに書き換えて
   per-face phin を独立な未知数として解くが、Dirichlet 側では
   対応する式の書き換え（per-face phin を未知数化）が未実装。
   書き換えなしで per-face Dirichlet DOF を振ると、同一点の phi が
   一意であるため BIE の行が線形従属になり行列が特異になる
  （実験で確認済み）。そのため、Neumann 面を持たない節点は棄却する。
--------------------------------------------------------------------------- */

inline bool hasSharpEdge(const auto* entity, double threshold = 60.0 * M_PI / 180.0) {
  return entity && entity->SharpQ(threshold);
}

inline bool detectMultipleNode(const auto* entity) {
  return hasSharpEdge(entity);
}

inline void setMultipleNode(auto* entity) {
  entity->isMultipleNode = detectMultipleNode(entity);

  // 棄却: Dirichlet のみの multiple node は現在の BIE 式では特異になる。
  // Neumann multiple node では式を書き換えて per-face phin を解くが、
  // Dirichlet 側の式書き換え（per-face phin 未知数化）は未実装。
  // 書き換えなしでは同一点の phi が一意なため行列が特異になる（実験確認済み）。
  if (entity->isMultipleNode && !hasAnyNeumannBoundaryState(entity))
    entity->isMultipleNode = false;
}
