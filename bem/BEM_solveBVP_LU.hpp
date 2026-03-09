#pragma once

// NOTE: This header is included inside `struct BEM_BVP` (see BEM_solveBVP.hpp).
// LU-specific (dense matrix) assembly.

void setIGIGn() {

  for (const auto &water : WATERS) {
    water->setGeometricPropertiesForce();
#pragma omp parallel for
    for (const auto &integ_f : water->getBoundaryFaces())
      integ_f->setIntegrationInfo();
  }

  this->IGIGn.assign(this->matrix_size, std::vector<std::array<double, 2>>(this->matrix_size, {0., 0.}));

  // TimeWatch timer;
  std::cout << "原点を節点にとり，方程式を作成．並列化" << std::endl;

  /*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

  #### 係数行列の作成

  数値シミュレーションでは，境界値問題を${\bf A}{\bf x}={\bf b}$のような線形連立方程式になるよう近似，変形し（離散化），${\bf x}$を求めることが多い．
  BEMでもBIEを離散化してこのような形にする．その際，境界条件に応じて，方程式（${\bf A}{\bf x}={\bf b}$の行）の右辺と左辺が入れ替える必要があるので注意する．
  これは，${\bf A}{\bf x}={\bf b}$の未知変数${\bf x}$と既知変数${\bf b}$がポテンシャル$\phi$か法線方向のポテンシャル$\phi_n$か，境界条件によって違うからである．
  プログラム上では，係数行列$\bf A$やベクトル$\bf b$を境界条件に応じて適切に作成すれば，求まる$\bf x$が適切なものになる．

  $\phi$の係数行列を$\mathbf{M}$，$\phi_n$の係数行列を$\mathbf{N}$，$\mathbf{\Phi}$を$\phi$のベクトル，$\mathbf{\Phi_n}$を$\phi_n$のベクトルとして，
  次のような連立一次方程式を得る．

  $$
  \mathbf{N} \mathbf{\Phi_n} = \mathbf{M} \mathbf{\Phi} \rightarrow {\bf A}{\bf x}={\bf b}
  $$

  このプログラムでは，$A$を`IGIGn`，$b$を`knowns`としている．

  このループでは，BIEの連立一次方程式の係数行列`IGIGn`を作成する作業を行なっている．
  `IGIGn`は，ある節点$i_\circ$（係数行列の行インデックス）に対する
  他の節点$j_\circ$（係数行列の列インデックス）の影響度合いのようなものである．
  その影響度合いは，他の節点$j_\circ$の所属する要素までの距離や向きによって決まることが離散化された式からわかる．

  | Variable | Description |
  |:--------:|:-----------:|
  | `origin` | 原点となる節点$i_\circ$ |
  | `integ_f` | Element $k_{\triangle}$ |
  | `t0, t1, ww` | Gaussian points and thier wieghts $\xi_0, \xi_1, w_0 w_1$ |
  | `p0, p1, p2` | Node of the element $k_{\triangle}$ |
  | `N012` | Shape function $\pmb{N}_j$ |
  | `IGIGn` | Coefficient matrices of the left and right sides |
  | `nr` | $\| \pmb{x} - \pmb{x}_{i\circ } \|$ |
  | `tmp` | $w_0 w_1 \frac{1 - \xi_0}{\| \pmb{x} - \pmb{x}_{i\circ } \|}$ |
  | `cross` | $\frac{\partial \pmb{x}}{\partial \xi_0} \times \frac{\partial \pmb{x}}{\partial \xi_1}$ |

  */

  for (const auto &water : WATERS) {
    water->setFacesVector();
    water->setPointsVector();
  }

  const T3Tdd shape_map_center = {{{0., 0.5} /*quad 4 -> linear 0*/, {0.5, 0.} /*quad 5 -> linear 1*/, {0.5, 0.5} /*quad 3 -> linear 2*/}};
  //! これとN3の内積を新なパラメタとして利用すると，(t0,t1)=(1,0)で{0., 0.5}に，(t0,t1)=(0,1)で{0.5, 0.}に，t0=t1=0で{0.5, 0.5}になる．
  const T3Tdd shape_map_l0_face = {{{0.5, 0.} /*0*/, {0., 0.5} /*1*/, {0., 0.} /*2*/}};
  const T3Tdd shape_map_l1_face = {{{0., 0.} /*0*/, {0.5, 0.} /*1*/, {0., 0.5} /*2*/}};
  const T3Tdd shape_map_l2_face = {{{0., 0.5} /*0*/, {0., 0.} /*1*/, {0.5, 0.} /*2*/}};

  constexpr std::array<double, 2> ZEROS2 = {0., 0.};
  constexpr std::array<double, 3> ZEROS3 = {0., 0., 0.};

  // // auto t0_t1_ww_N012_HIGHRESOLUTION = t0_t1_ww_N012_LOWRESOLUTION;

  // for (int i = 0; const auto &[t0, t1, ww] : __array_GW10xGW10__) {
  //    auto t0t1t2 = ModTriShape<3>(t0, t1);
  //    t0_t1_ww_N012_HIGHRESOLUTION.push_back({t0, t1, ww, t0t1t2});
  // }

  /*
  ## 積分の効率化

  `linear_triangle_integration_info`と`pseudo_quadratic_triangle_integration_info`は，
  積分点の位置と重みを事前に計算しておくことで，積分の効率化を図るためのものである．
  `linear_triangle_integration_info`と`pseudo_quadratic_triangle_integration_info`の引数の詳細

  ### `linear_triangle_integration_info`
  | 引数名  | 説明                                                              |
  |---------|-----------------------------------------------------------------|
  | `Tdd`   | 2D パラメータ {[0,1], [0,1]}：積分変数                             |
  | `double`| ガウス重み（積分重み）                                             |
  | `Tddd`  | 3D パラメータ {xi0=[0,1], xi1=[0,1-xi0]}：2Dパラメータと関連        |
  | `Tddd`  | 3D 位置ベクトル（{xi0, xi1, xi2}を使用）                           |
  | `Tddd`  | 外積（dX/dxi0 × dX/dxi1）                                        |
  | `double`| 外積のノルム                                                     |

  ### `pseudo_quadratic_triangle_integration_info`
  | 引数名              | 説明                                                              |
  |---------------------|-----------------------------------------------------------------|
  | `Tdd`               | 2D パラメータ {[0,1], [0,1]}：積分変数                             |
  | `double`            | ガウス重み（積分重み）                                             |
  | `Tddd`              | 3D パラメータ {xi0=[0,1], xi1=[0,1-xi0]}：2Dパラメータと関連        |
  | `std::array<T6d, 4>`| 二次要素の形状関数                                                |
  | `Tddd`              | 3D 位置ベクトル（{xi0, xi1, xi2}を使用）                           |
  | `Tddd`              | 外積（dX/dxi0 × dX/dxi1）                                        |
  | `double`            | 外積のノルム                                                     |
  */

  int count_pseudo_quadratic_element = 0, count_linear_element = 0, total = 0;
  for (const auto &water : WATERS)
    for (const auto &integ_f : water->getBoundaryFaces()) {
      if (integ_f->isLinearElement)
        count_linear_element++;
      else if (integ_f->isPseudoQuadraticElement)
        count_pseudo_quadratic_element++;
      total++;
    }

  std::cout << "線形要素の面の数：" << count_linear_element << " persecent: " << 100. * count_linear_element / total << std::endl;
  std::cout << "擬似二次要素の面の数：" << count_pseudo_quadratic_element << " persecent: " << 100. * count_pseudo_quadratic_element / total << std::endl;

  if (use_pseudo_quadratic_element)
    std::cout << "擬似二次要素を使ってBIEを離散化" << std::endl;
  else
    std::cout << "線形要素を使ってBIEを離散化" << std::endl;

  for (const auto water : WATERS) {
    double scale = water->getScale();
    auto surfacePoints = water->getBoundaryPoints();
    auto surfaces = water->getBoundaryFaces();
#pragma omp parallel for
    for (const auto &origin : surfacePoints) {
      //@ this loop is for the multiple nodes
      for (const auto &[_, index] : origin->f2Index) {
        double origin_ign_rigid_mode = 0.;
        auto &IGIGn_Row = IGIGn[index];
        double nr, ig, ign, ww_nr, tmp;
        Tdd ig_ign0, ig_ign1, ig_ign2;
        Tddd R;
        networkPoint *closest_p_to_origin = nullptr;
        std::vector<std::tuple<networkPoint *, networkFace *, double, double>> key_ig_ign;

        //@ for all water faces
        for (const auto &water : WATERS) {
          // b@ integrate over all faces
          for (const auto &integ_f : surfaces) {
            auto [p0, p1, p2] = integ_f->getPoints(origin);
            closest_p_to_origin = p0;

            /*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

            WARNING: この`std::vector<std::tuple<networkPoint *, networkFace *, double, double>> key_ig_ign`の`networkFace`は，どの面側から節点を呼び出すかを決めていて，高次補間の場合，積分面と一致しない場合がある．

            1. fill key_ig_ign
            2. fill IGIGn_Row

            */
            auto dist = Norm((p0->X + p1->X + p2->X) / 3. - origin->X);
            int how_far = 0;
            if (dist < scale / 25.)
              how_far = 1;

            if (integ_f->isLinearElement) {
              ig_ign0 = ig_ign1 = ig_ign2 = {0., 0.};
              if (p0 == origin) {
                ign = 0.;
                for (const auto &[t0t1, ww, shape3, X, cross, J_det] : integ_f->map_Point_LinearIntegrationInfo_vector[1].at(closest_p_to_origin)) {
                  ig = J_det * (ww / (nr = Norm(R = (X - origin->X))));
                  std::get<0>(ig_ign0) += ig * std::get<0>(shape3);
                  std::get<0>(ig_ign1) += ig * std::get<1>(shape3);
                  std::get<0>(ig_ign2) += ig * std::get<2>(shape3);
                }
              } else {
                for (const auto &[t0t1, ww, shape3, X, cross, J_det] : integ_f->map_Point_LinearIntegrationInfo_vector[how_far].at(closest_p_to_origin)) {
                  ig = J_det * (ww_nr = ww / (nr = Norm(R = (X - origin->X))));
                  ign = Dot(R, cross) * ww_nr / (nr * nr);
                  std::get<0>(ig_ign0) += ig * std::get<0>(shape3);
                  std::get<1>(ig_ign0) -= ign * std::get<0>(shape3);
                  std::get<0>(ig_ign1) += ig * std::get<1>(shape3);
                  std::get<1>(ig_ign1) -= ign * std::get<1>(shape3);
                  std::get<0>(ig_ign2) += ig * std::get<2>(shape3);
                  std::get<1>(ig_ign2) -= ign * std::get<2>(shape3);
                  origin_ign_rigid_mode += ign;
                }
              }
              IGIGn_Row[pf2Index(p0, integ_f)] += ig_ign0;
              IGIGn_Row[pf2Index(p1, integ_f)] += ig_ign1;
              IGIGn_Row[pf2Index(p2, integ_f)] += ig_ign2;
            } else if (integ_f->isPseudoQuadraticElement) {
              key_ig_ign = integ_f->map_Point_BEM_IGIGn_info_init.at(closest_p_to_origin);
              for (const auto &[t0t1, ww, Nc_N0_N1_N2, X, cross, J_det] : integ_f->map_Point_PseudoQuadraticIntegrationInfo_vector[how_far].at(closest_p_to_origin)) {
                ig = J_det * (ww_nr = ww / (nr = Norm(R = (X - origin->X))));
                ign = Dot(R, cross) * ww_nr / (nr * nr);
                for (auto i = 0; i < 6; ++i) {
                  tmp = std::get<0>(Nc_N0_N1_N2)[i];
                  std::get<2>(key_ig_ign[i]) += ig * tmp;
                  std::get<3>(key_ig_ign[i]) -= ign * tmp;
                  if (std::get<0>(key_ig_ign[i]) != origin)
                    origin_ign_rigid_mode += ign * tmp;

                  tmp = std::get<1>(Nc_N0_N1_N2)[i];
                  std::get<2>(key_ig_ign[i + 6]) += ig * tmp;
                  std::get<3>(key_ig_ign[i + 6]) -= ign * tmp;
                  if (std::get<0>(key_ig_ign[i + 6]) != origin)
                    origin_ign_rigid_mode += ign * tmp;

                  tmp = std::get<2>(Nc_N0_N1_N2)[i];
                  std::get<2>(key_ig_ign[i + 12]) += ig * tmp;
                  std::get<3>(key_ig_ign[i + 12]) -= ign * tmp;
                  if (std::get<0>(key_ig_ign[i + 12]) != origin)
                    origin_ign_rigid_mode += ign * tmp;

                  tmp = std::get<3>(Nc_N0_N1_N2)[i];
                  std::get<2>(key_ig_ign[i + 18]) += ig * tmp;
                  std::get<3>(key_ig_ign[i + 18]) -= ign * tmp;
                  if (std::get<0>(key_ig_ign[i + 18]) != origin)
                    origin_ign_rigid_mode += ign * tmp;
                }
              }

              for (const auto &[p, integ_f, ig, ign] : key_ig_ign)
                IGIGn_Row[pf2Index(p, integ_f)] += Tdd{ig, ign}; // この面に関する積分において，φまたはφnの寄与
            } else if (integ_f->isTrueQuadraticElement) {
              // True quadratic integration (vertex origin)
              auto [fp0, fl0, fp1, fl1, fp2, fl2] = integ_f->PLPLPL;
              constexpr std::array<bool, 3> all_true{true, true, true};
              std::array<Tddd, 6> X6 = {fp0->X, fp1->X, fp2->X, fl0->X_mid, fl1->X_mid, fl2->X_mid};
              auto dof_indices = getQuadDOFIndices(integ_f);

              std::array<std::array<double, 2>, 6> WGN_WGnN{};

              auto redistribute_and_flush = [&]() {
                // Redistribute inactive midpoint contributions to endpoint vertices
                for (int k = 3; k < 6; ++k) {
                  if (dof_indices[k] < 0) {
                    int v1 = k - 3, v2 = (k - 3 + 1) % 3;
                    WGN_WGnN[v1][0] += 0.5 * WGN_WGnN[k][0];
                    WGN_WGnN[v1][1] += 0.5 * WGN_WGnN[k][1];
                    WGN_WGnN[v2][0] += 0.5 * WGN_WGnN[k][0];
                    WGN_WGnN[v2][1] += 0.5 * WGN_WGnN[k][1];
                    WGN_WGnN[k] = {0., 0.};
                  }
                }
                for (int k = 0; k < 6; ++k)
                  if (dof_indices[k] >= 0) {
                    IGIGn_Row[dof_indices[k]] += Tdd{WGN_WGnN[k][0], WGN_WGnN[k][1]};
                    if (dof_indices[k] != index)
                      origin_ign_rigid_mode += WGN_WGnN[k][1];
                  }
              };

              // Check if origin is a vertex of this face
              int vertex_local = -1;
              if (fp0 == origin) vertex_local = 0;
              else if (fp1 == origin) vertex_local = 1;
              else if (fp2 == origin) vertex_local = 2;

              if (vertex_local >= 0) {
                // --- Vertex Duffy ---
                std::array<int, 3> vperm = {0, 1, 2};
                if (vertex_local == 1) vperm = {1, 2, 0};
                else if (vertex_local == 2) vperm = {2, 0, 1};

                for (const auto &[t0, t1, ww] : __array_GW5xGW5__) {
                  auto bary_loc = ModTriShape<3>(t0, t1);
                  std::array<double, 3> bary{};
                  bary[vperm[0]] = bary_loc[0];
                  bary[vperm[1]] = bary_loc[1];
                  bary[vperm[2]] = bary_loc[2];

                  auto N6_geo = TriShape<6>(bary[0], bary[1], all_true);
                  auto dN_dt0 = D_TriShape<6, 1, 0>(bary[0], bary[1], all_true);
                  auto dN_dt1 = D_TriShape<6, 0, 1>(bary[0], bary[1], all_true);
                  auto X_q = Dot(N6_geo, X6);
                  auto cross_q = Cross(Dot(dN_dt0, X6), Dot(dN_dt1, X6));
                  auto N6 = integ_f->trueQuadN6(bary[0], bary[1]);

                  double nr_v = Norm(R = (X_q - origin->X));
                  double tmp_v = ww * (1. - t0) / nr_v;
                  double WG = Norm(cross_q) * tmp_v;
                  double WGn = -Dot(R / nr_v, cross_q) * tmp_v / nr_v;
                  for (int k = 0; k < 6; ++k) {
                    WGN_WGnN[k][0] += WG * N6[k];
                    WGN_WGnN[k][1] += (k == vertex_local) ? 0. : WGn * N6[k];
                  }
                }
                redistribute_and_flush();
              } else {
                // --- Non-adjacent (vertex origin not on this face) ---
                const bool is_near = (dist < scale / 25.);
                if (is_near) {
                  for (const auto &[t0, t1, ww] : __array_GW5xGW5__) {
                    auto bary_loc = ModTriShape<3>(t0, t1);
                    auto N6_geo = TriShape<6>(bary_loc[0], bary_loc[1], all_true);
                    auto dN_dt0 = D_TriShape<6, 1, 0>(bary_loc[0], bary_loc[1], all_true);
                    auto dN_dt1 = D_TriShape<6, 0, 1>(bary_loc[0], bary_loc[1], all_true);
                    auto X_q = Dot(N6_geo, X6);
                    auto cross_q = Cross(Dot(dN_dt0, X6), Dot(dN_dt1, X6));
                    auto N6 = integ_f->trueQuadN6(bary_loc[0], bary_loc[1]);
                    double weight = ww * (1. - t0);

                    double nr_v = Norm(R = (X_q - origin->X));
                    double nr_inv = 1. / nr_v;
                    double WG = Norm(cross_q) * weight * nr_inv;
                    double WGn = -Dot(R * nr_inv, cross_q) * weight * nr_inv * nr_inv;
                    for (int k = 0; k < 6; ++k) {
                      WGN_WGnN[k][0] += WG * N6[k];
                      WGN_WGnN[k][1] += WGn * N6[k];
                    }
                    origin_ign_rigid_mode += WGn;
                  }
                } else {
                  for (const auto &[b0, b1, b2, ww] : getDunavantP2MRule(6)) {
                    auto N6_geo = TriShape<6>(b0, b1, all_true);
                    auto dN_dt0 = D_TriShape<6, 1, 0>(b0, b1, all_true);
                    auto dN_dt1 = D_TriShape<6, 0, 1>(b0, b1, all_true);
                    auto X_q = Dot(N6_geo, X6);
                    auto cross_q = Cross(Dot(dN_dt0, X6), Dot(dN_dt1, X6));
                    auto N6 = integ_f->trueQuadN6(b0, b1);

                    double nr_v = Norm(R = (X_q - origin->X));
                    double nr_inv = 1. / nr_v;
                    double WG = Norm(cross_q) * ww * nr_inv;
                    double WGn = -Dot(R * nr_inv, cross_q) * ww * nr_inv * nr_inv;
                    for (int k = 0; k < 6; ++k) {
                      WGN_WGnN[k][0] += WG * N6[k];
                      WGN_WGnN[k][1] += WGn * N6[k];
                    }
                    origin_ign_rigid_mode += WGn;
                  }
                }
                redistribute_and_flush();
              }
            }
          }
        }

        /* -------------------------------------------------------------------------- */
        /*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

        ### リジッドモードテクニック（係数行列の対角成分の計算）

        BIEの対角成分の計算で注意が必要なのは，原点$i_\circ$の頂点の立体角と，係数の特異性である．

        * 係数行列の対角成分には，立体角$\alpha$が含まれており，この計算は面倒である．
        * 係数の計算には，$\frac{{\mathbf{x}_{k_\vartriangle}({\pmb{\xi}}) - \mathbf{x}_{i_\circ}}}{{\| \mathbf{x}_{k_\vartriangle}({\pmb{\xi}}) - \mathbf{x}_{i_\circ} \|}^3}$が含まれており，分母が0付近で強い特異性を持つ．

        そこで，素直に幾何学的な観点から立体角を計算するのではなく，BIEの式を使って積分で計算する方法がある．BIEの式に，$\phi=1$を代入すると，$\phi_n$が消える．結局，対角成分，つまり，原点$i_\circ$を頂点上の変数に掛かる係数は，次のようになる．

        $$
        \sum\limits_{k_\vartriangle} 2 A_{k_\vartriangle} \, \mathbf{n}_{k_\vartriangle} \cdot \sum\limits_{\xi_1, w_1} \sum\limits_{\xi_0, w_0} \left( w_0 w_1 \left( \sum\limits_{j=0}^2 \bar\delta_{(k_\vartriangle, j),i_\circ} N_j({\pmb{\xi}}) \right) \frac{{\mathbf{x}_{k_\vartriangle}({\pmb{\xi}}) - \mathbf{x}_{i_\circ}}}{{\| \mathbf{x}_{k_\vartriangle}({\pmb{\xi}}) - \mathbf{x}_{i_\circ} \|}^3}(1 - \xi_0)\right)
        $$

        $\bar\delta_{(k_\vartriangle, j),i_\circ}$は，$k_\vartriangle$の$j$番目の頂点が$i_\circ$である場合に0，それ以外は1となる関数である．

        数値計算上は，$\delta_{(k_\vartriangle, j),i_\circ}$がゼロの場合は，そもそも係数をインクリメントせず，スキップする．
        これはリジッドモードテクニックと呼ばれていて，分子が小さくなる特異的な計算を省き，立体角の計算もまとめて対角成分を計算することができる方法である．

        ただし，線形要素の場合，原点$i_\circ$を頂点とする三角形$k_{\vartriangle}$に対する計算，${\bf n}_{k_\vartriangle}\cdot ({{\bf x}_{k_\vartriangle}}(\pmb{\xi})-{{\bf x}_{i_\circ}})=0$となるため，和をとる必要はない．
        よって，そもそも線形要素の場合は，特異的な計算は含まれない．

        */

        // #if defined(use_rigid_mode)
        std::get<1>(IGIGn_Row[index]) = origin_ign_rigid_mode;
        // #else
        //           std::get<1>(IGIGn_Row[index]) += origin->getSolidAngle();
        // #endif
      }
    }

    // Midpoint origin loop for true quadratic elements
    if (use_true_quadratic_element) {
      auto midpoint_lines = water->getBoundaryLines();
#pragma omp parallel for
      for (const auto &origin_line : midpoint_lines) {
        for (const auto &[_, index] : origin_line->f2Index) {
          double origin_ign_rigid_mode = 0.;
          auto &IGIGn_Row = IGIGn[index];
          Tddd R;
          const Tddd origin_pos = origin_line->Xtarget;

          for (const auto &water : WATERS) {
            for (const auto &integ_f : surfaces) {
              if (!integ_f->isTrueQuadraticElement)
                continue;

              auto [fp0, fl0, fp1, fl1, fp2, fl2] = integ_f->PLPLPL;
              constexpr std::array<bool, 3> all_true{true, true, true};
              std::array<Tddd, 6> X6 = {fp0->X, fp1->X, fp2->X, fl0->X_mid, fl1->X_mid, fl2->X_mid};
              auto dof_indices = getQuadDOFIndices(integ_f);
              auto dist = Norm(integ_f->centroid - origin_pos);

              std::array<std::array<double, 2>, 6> WGN_WGnN{};

              auto redistribute_and_flush = [&]() {
                for (int k = 3; k < 6; ++k) {
                  if (dof_indices[k] < 0) {
                    int v1 = k - 3, v2 = (k - 3 + 1) % 3;
                    WGN_WGnN[v1][0] += 0.5 * WGN_WGnN[k][0];
                    WGN_WGnN[v1][1] += 0.5 * WGN_WGnN[k][1];
                    WGN_WGnN[v2][0] += 0.5 * WGN_WGnN[k][0];
                    WGN_WGnN[v2][1] += 0.5 * WGN_WGnN[k][1];
                    WGN_WGnN[k] = {0., 0.};
                  }
                }
                for (int k = 0; k < 6; ++k)
                  if (dof_indices[k] >= 0) {
                    IGIGn_Row[dof_indices[k]] += Tdd{WGN_WGnN[k][0], WGN_WGnN[k][1]};
                    if (dof_indices[k] != index)
                      origin_ign_rigid_mode += WGN_WGnN[k][1];
                  }
              };

              // Check if origin_line is an edge midpoint of this face
              int on_edge = -1;
              if (fl0 == origin_line) on_edge = 0;
              else if (fl1 == origin_line) on_edge = 1;
              else if (fl2 == origin_line) on_edge = 2;

              if (on_edge >= 0) {
                // --- Edge midpoint Duffy (2 sub-triangles) ---
                int va = on_edge, vb = (on_edge + 1) % 3, vc = (on_edge + 2) % 3;
                const int mid_local = 3 + on_edge;

                for (const auto &[t0, t1, ww] : __array_GW5xGW5__) {
                  double omt0 = 1.0 - t0;
                  double s1 = t1 * omt0;
                  double s2 = omt0 * (1.0 - t1);

                  // Sub-triangle 1: (M, vb, vc)
                  {
                    std::array<double, 3> bary{};
                    bary[va] = 0.5 * t0;
                    bary[vb] = 0.5 * t0 + s1;
                    bary[vc] = s2;
                    auto N6_geo = TriShape<6>(bary[0], bary[1], all_true);
                    auto dN_dt0 = D_TriShape<6, 1, 0>(bary[0], bary[1], all_true);
                    auto dN_dt1 = D_TriShape<6, 0, 1>(bary[0], bary[1], all_true);
                    auto X_q = Dot(N6_geo, X6);
                    auto cross_q = Cross(Dot(dN_dt0, X6), Dot(dN_dt1, X6));
                    auto N6 = integ_f->trueQuadN6(bary[0], bary[1]);

                    double nr_v = Norm(R = (X_q - origin_pos));
                    double tmp_v = ww * omt0 * 0.5 / nr_v;
                    double WG = Norm(cross_q) * tmp_v;
                    double WGn = -Dot(R / nr_v, cross_q) * tmp_v / nr_v;
                    for (int k = 0; k < 6; ++k) {
                      WGN_WGnN[k][0] += WG * N6[k];
                      WGN_WGnN[k][1] += (k == mid_local) ? 0. : WGn * N6[k];
                    }
                  }

                  // Sub-triangle 2: (M, vc, va)
                  {
                    std::array<double, 3> bary{};
                    bary[va] = 0.5 * t0 + s2;
                    bary[vb] = 0.5 * t0;
                    bary[vc] = s1;
                    auto N6_geo = TriShape<6>(bary[0], bary[1], all_true);
                    auto dN_dt0 = D_TriShape<6, 1, 0>(bary[0], bary[1], all_true);
                    auto dN_dt1 = D_TriShape<6, 0, 1>(bary[0], bary[1], all_true);
                    auto X_q = Dot(N6_geo, X6);
                    auto cross_q = Cross(Dot(dN_dt0, X6), Dot(dN_dt1, X6));
                    auto N6 = integ_f->trueQuadN6(bary[0], bary[1]);

                    double nr_v = Norm(R = (X_q - origin_pos));
                    double tmp_v = ww * omt0 * 0.5 / nr_v;
                    double WG = Norm(cross_q) * tmp_v;
                    double WGn = -Dot(R / nr_v, cross_q) * tmp_v / nr_v;
                    for (int k = 0; k < 6; ++k) {
                      WGN_WGnN[k][0] += WG * N6[k];
                      WGN_WGnN[k][1] += (k == mid_local) ? 0. : WGn * N6[k];
                    }
                  }
                }
                redistribute_and_flush();
              } else {
                // --- Non-adjacent (midpoint origin not on this face) ---
                const bool is_near = (dist < scale / 25.);
                if (is_near) {
                  for (const auto &[t0, t1, ww] : __array_GW5xGW5__) {
                    auto bary_loc = ModTriShape<3>(t0, t1);
                    auto N6_geo = TriShape<6>(bary_loc[0], bary_loc[1], all_true);
                    auto dN_dt0 = D_TriShape<6, 1, 0>(bary_loc[0], bary_loc[1], all_true);
                    auto dN_dt1 = D_TriShape<6, 0, 1>(bary_loc[0], bary_loc[1], all_true);
                    auto X_q = Dot(N6_geo, X6);
                    auto cross_q = Cross(Dot(dN_dt0, X6), Dot(dN_dt1, X6));
                    auto N6 = integ_f->trueQuadN6(bary_loc[0], bary_loc[1]);
                    double weight = ww * (1. - t0);

                    double nr_v = Norm(R = (X_q - origin_pos));
                    double nr_inv = 1. / nr_v;
                    double WG = Norm(cross_q) * weight * nr_inv;
                    double WGn = -Dot(R * nr_inv, cross_q) * weight * nr_inv * nr_inv;
                    for (int k = 0; k < 6; ++k) {
                      WGN_WGnN[k][0] += WG * N6[k];
                      WGN_WGnN[k][1] += WGn * N6[k];
                    }
                    origin_ign_rigid_mode += WGn;
                  }
                } else {
                  for (const auto &[b0, b1, b2, ww] : getDunavantP2MRule(6)) {
                    auto N6_geo = TriShape<6>(b0, b1, all_true);
                    auto dN_dt0 = D_TriShape<6, 1, 0>(b0, b1, all_true);
                    auto dN_dt1 = D_TriShape<6, 0, 1>(b0, b1, all_true);
                    auto X_q = Dot(N6_geo, X6);
                    auto cross_q = Cross(Dot(dN_dt0, X6), Dot(dN_dt1, X6));
                    auto N6 = integ_f->trueQuadN6(b0, b1);

                    double nr_v = Norm(R = (X_q - origin_pos));
                    double nr_inv = 1. / nr_v;
                    double WG = Norm(cross_q) * ww * nr_inv;
                    double WGn = -Dot(R * nr_inv, cross_q) * ww * nr_inv * nr_inv;
                    for (int k = 0; k < 6; ++k) {
                      WGN_WGnN[k][0] += WG * N6[k];
                      WGN_WGnN[k][1] += WGn * N6[k];
                    }
                    origin_ign_rigid_mode += WGn;
                  }
                }
                redistribute_and_flush();
              }
            }
          }

          // Rigid mode for midpoint origin
          std::get<1>(IGIGn_Row[index]) = origin_ign_rigid_mode;
        }
      }
    }
  }
  // std::cout << Green << "離散化にかかった時間" << timer() << colorReset << std::endl;
};

void generateBIEMatrix() {
  std::cout << "generateBIEMatrix()" << std::endl;

  /*DOC_EXTRACT 0_2_BOUNDARY_VALUE_PROBLEM

  ### 左辺と右辺の入れ替え

  係数行列`IGIGn`は，左辺の$I_G \phi_n$，右辺の$I_{G_n}\phi$の係数行列を表している．

  $$
  (I_G)_{i_\circ,j_\circ} (\phi_n)_{j_\circ} = (I_{Gn})_{i_\circ,j_\circ}  \phi_{j_\circ}
  $$

  境界条件に応じて，未知変数は$\phi,\phi_n$のどちらかに決まる．
  未知変数が$\phi$の場合（Dirichlet境界条件の場合），
  係数行列`IGIGn`中で対応する列を符号変えて入れ替えることで移項したことになる．

  #### ２種類の多重節点

  1. Dirichlet面上であり，かつNeumann面上である多重節点
  2. Dirichlet面上ではなく，完全にNeumann面上にあるが，法線ベクトルが大きく異なる節点

  1の多重節点の場合，BIEの連立一次方程式の係数行列の行を，Dirchlet面上の$\phi$とNeumann面上の$\phi$の値が一致する，という式に変更する．
  2の場合は，特に変更しない．BIEを解くことで，それぞれの面に対して，$\phi$が得られるが，それらの平均値，または重み付け平均値を$\phi$として採用する．

  */

  //^ ---------------------------- setIGIGnで計算した結果を基に，係数行列を作成する．---------------------------
  double max_value = 0;
  for (auto i = 0; i < IGIGn.size(); ++i) {
    if (max_value < std::abs(std::get<0>(IGIGn[i][i])))
      max_value = std::abs(std::get<0>(IGIGn[i][i]));
  }

  this->mat_kn.assign(this->matrix_size, V_d(this->matrix_size, 0.));
  this->mat_ukn.assign(this->matrix_size, V_d(this->matrix_size, 0.));

  // Helper lambda: fill columns for a given row index i from IGIGn
  auto fill_columns = [&](int i, const auto &surfacePoints, const auto &water) {
    auto &IGIGn_i = IGIGn[i];
    // Vertex DOF columns
    for (const auto &x : surfacePoints)
      for (const auto &[x_face, j] : x->f2Index) {
        mat_ukn[i][j] = IGIGn_i[j][0];
        mat_kn[i][j] = IGIGn_i[j][1];
        if (isNeumannID_BEM(x, x_face)) {
          std::swap(mat_ukn[i][j], mat_kn[i][j]);
          mat_ukn[i][j] *= -1;
          mat_kn[i][j] *= -1;
        }
      }
    // Midpoint DOF columns (true quadratic only)
    if (use_true_quadratic_element)
      for (auto *x_line : water->getBoundaryLines())
        for (const auto &[x_face, j] : x_line->f2Index) {
          mat_ukn[i][j] = IGIGn_i[j][0];
          mat_kn[i][j] = IGIGn_i[j][1];
          if (isNeumannID_BEM(x_line, x_face)) {
            std::swap(mat_ukn[i][j], mat_kn[i][j]);
            mat_ukn[i][j] *= -1;
            mat_kn[i][j] *= -1;
          }
        }
  };

  for (const auto water : WATERS) {
    auto surfacePoints = water->getBoundaryPoints();
    // Vertex DOF rows
#pragma omp parallel for
    for (const auto &a : surfacePoints)
      for (const auto &[a_face, i] : a->f2Index) {
        fill_columns(i, surfacePoints, water);

        //^ 多重節点の場合，Neumann面上のphiの値は，同じ場所のDirichlet面のphiと一致する
        if (a->CORNER && isNeumannID_BEM(a, a_face) /*行の変更*/) {
          std::ranges::fill(mat_ukn[i], 0.);
          std::ranges::fill(mat_kn[i], 0.);
          mat_ukn[i][i] = max_value;
          mat_kn[i][pf2Index(a, nullptr)] = max_value;
        }
      }

    // Midpoint DOF rows (true quadratic only)
    if (use_true_quadratic_element) {
#pragma omp parallel for
      for (auto *a_line : water->getBoundaryLines())
        for (const auto &[a_face, i] : a_line->f2Index) {
          fill_columns(i, surfacePoints, water);

          // CORNER midpoint: Neumann face phi = Dirichlet face phi
          if (a_line->CORNER && isNeumannID_BEM(a_line, a_face)) {
            std::ranges::fill(mat_ukn[i], 0.);
            std::ranges::fill(mat_kn[i], 0.);
            mat_ukn[i][i] = max_value;
            mat_kn[i][lf2Index(a_line, nullptr)] = max_value;
          }
        }
    }
  }
  knowns = getVectorFromBoundary(WATERS, this->matrix_size, [](const networkPoint *p, networkFace *f) { return p->phiOnFace.at(f); }, [](const networkPoint *p, networkFace *f) { return p->phinOnFace.at(f); });
  b_RHS = Dot(mat_kn, knowns);
};

void solveLU(TimeWatch &watch, double &time_setup, double &time_solve) {
  std::cout << "Using LU solver..." << std::endl;
  watch();
  setIGIGn();
  // std::cout << "2つの係数行列の情報を持つ　P_P_IGIGn　を境界条件に応じて入れ替える（移項）:" << std::endl;
  generateBIEMatrix();
  std::cout << "BIE setup (dense matrix build) time: " << Green << (time_setup = watch()[0]) << colorReset << std::endl;
  ans.resize(knowns.size(), 0.);
  //! A.x = b
  if (this->lu != nullptr) {
    this->lu->init(mat_ukn);
    this->lu->solve(b_RHS /*既知のベクトル（右辺）*/, ans /*解*/);
  } else
    this->lu = new lapack_lu(mat_ukn /*未知の行列係数（左辺）*/, ans /*解*/, b_RHS /*既知のベクトル（右辺）*/);
  const double time_lapack_lu = watch()[0];
  std::cout << "LU分解を使った解法に要した時間: " << Green << time_lapack_lu << colorReset << std::endl;
  time_solve = time_lapack_lu;
  std::cout << colorReset << "update p->phiphin and p->phinOnFace for Dirichlet boundary" << colorReset << std::endl;
  isSolutionFinite(ans);
  storePhiPhin(WATERS, ans);
}
