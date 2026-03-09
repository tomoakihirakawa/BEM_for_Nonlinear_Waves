#ifndef calcDerivatives_H
#define calcDerivatives_H
#pragma once
////////////////////////////////////////////////////

template <typename T, typename U> std::unordered_map<T, U> map2unordered_map(const std::map<T, U> map_IN) {
  std::unordered_map<T, U> ret;
  ret.reserve(map_IN.size());
  for (const auto &[a, b] : map_IN)
    ret[a] = b;
  return ret;
};
template <typename T, typename U> std::unordered_map<T *, U> map2unordered_map(const std::map<T *, U> map_IN) {
  std::unordered_map<T *, U> ret;
  ret.reserve(map_IN.size());
  for (const auto &[a, b] : map_IN)
    ret[a] = b;
  return ret;
};
////////////////////////////////////////////////////
class calcDerivatives {
  using map_P_d = std::map<netP *, double>;
  using map_P_Vd = std::map<netP *, V_d>;
  using map_F_P_Vd = std::map<netF *, map_P_Vd>;
  using map_P_P_Vd = std::map<netP *, map_P_Vd>;
  using map_P_F_P_Vd = std::map<netP *, map_F_P_Vd>;

  using uo_map_P_d = std::unordered_map<netP *, double>;
  using uo_map_P_Vd = std::unordered_map<netP *, V_d>;
  // using uo_map_F_P_Vd = std::unordered_map<netF *, map_P_Vd>;
  // using uo_map_P_P_Vd = std::unordered_map<netP *, map_P_Vd>;
  // using uo_map_P_F_P_Vd = std::unordered_map<netP *, map_F_P_Vd>;

public:
  map_P_Vd P_phiphin_InnerOuterCornerP; // InnerP, OuterP, CornerP
  map_P_Vd P_phiphin_InnerOuterP;       // InnerP, OuterP
  map_P_Vd P_phiphin_BIE;               // InnerP, CornerP: P_phiphin_InnerCornerPに等しい
  map_P_Vd P_nablaPhi;
  map_P_d P_DphiDt;
  map_P_P_Vd P_P_IGIGn; // 途中でreverseしているので注意
  VV_d matOfKnowns, matOfUnknowns;
  V_netPp pOfknowns;
  ludcmp *lu;
  SVD *svd;
  std::string solver_name;

  bool isWellConditioned(const netPp p) {
    for (const auto &l : p->getLines()) {
      // auto ps = l->getPoints();
      // if (!AllTrue(ps))
      // {
      //     std::stringstream ss;
      //     ss << "ps = " << ps;
      //     throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str()));
      // }
      // if (ps[0] != p && ps[1] != p)
      // { //l is not connected to p
      //     geometry::Point_Line PL(p->getX(), ps[0]->getX(), ps[1]->getX());
      //     if (PL.d2line_segment < 1E-3)
      //         return false;
      // }

      ///////////////////////////
      if (l->length() < 1E-3)
        return false;
    }
    return true;
  };

  ~calcDerivatives() {
    if (this->lu != nullptr)
      delete (this->lu);
    if (this->svd != nullptr)
      delete (this->lu);
  };
  calcDerivatives(map_P_Vd &given_P_phiphin_InnerOuterP /*正確なphiphinを与える*/ /*BIEで用いないものも含まれている*/, const BEM::CompGrid &cpg) : P_phiphin_InnerOuterP({}), P_phiphin_InnerOuterCornerP({}), P_phiphin_BIE({}), P_nablaPhi({}), P_DphiDt({}), matOfKnowns({}), matOfUnknowns({}), pOfknowns({}), lu(nullptr), svd(nullptr), knowns_DN({}), solver_name("") {
    auto Faces = cpg.well_Faces;
    auto InnerP = cpg.InnerP;
    auto CornerP = cpg.CornerP;
    auto OuterP = cpg.OuterP;

    ////////////////////////

    /*与えられたphi, phinがこのBIEを解くために用いられるべき計算領域の変数かはわかっていない*/
    //-------------- InnerP & CornerP -> P_phiphin_BIE -----------------
    Print("`given_P_phiphin_InnerOuterP`内で計算領域内にある点：`InnerP`を抽出．これとCornerPが計算領域であり，BIEでD->N,N->Dを計算する．", Green);
    Print("(i) InnerPでも新たな点は，phiとphinは計算されていないので補間する", Magenta);

    auto nets = takeNetworks(TakeFirst(given_P_phiphin_InnerOuterP));
    std::cout << "nets = " << getNames(nets) << std::endl;
    Print(nets, Green);

    for (const auto &p : InnerP) {
      if (given_P_phiphin_InnerOuterP.find(p) != given_P_phiphin_InnerOuterP.end())
        this->P_phiphin_BIE[p] = given_P_phiphin_InnerOuterP[p];
      else {
        if (p->isD()) {
          mk_vtu("./vtu/new_point.vtu", {{p}});
          Print(getNames({p}));
          throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "Dirichlet条件の水粒子が内部に新たに入ってくることはありえない");
          // this->P_phiphin_BIE[p] = {BEM::InterpolationRBF_phi(p, given_P_phiphin_InnerOuterP, netsampler.net4phi(p)), 0.};
          // this->P_phiphin_BIE[p] = {BEM::InterpolationRBF_phi(p, given_P_phiphin_InnerOuterP, getReferenceNetwork(p, TakeFirst(given_P_phiphin_InnerOuterP))), 0.};
          this->P_phiphin_BIE[p] = {BEM::InterpolationRBF_phi(p, given_P_phiphin_InnerOuterP, nets), 0.};
        } else if (p->isN())
          this->P_phiphin_BIE[p] = {0., 0.};
      }
    }

    Print("(0) `P_phiphin_BIE`に角点を加える", Green);
    // for (const auto &p : CornerP)
    //     this->P_phiphin_BIE[p] = {0., 0.};

    V_netPp CornerP_good = {};
    V_netPp CornerP_bad = {};

    for (const auto &p : CornerP) {
      if (isWellConditioned(p)) {
        this->P_phiphin_BIE[p] = {0., 0.};
        CornerP_good.emplace_back(p);
      } else
        CornerP_bad.emplace_back(p);
    }

    BEM::networkSampler netsamp(Join(takeNetworks(InnerP), takeNetworks(CornerP_good)));

    //============================================
    //================= IG IGn ===================
    //============================================
    /// calculate IG and IGn that satisfy
    /// IG *  phi_n  = (IGn - c * delta)  * phi

    Print("(1-1) P_P_IGIGnの計算", Green);

    V_netPp InnerCornerP = Join(InnerP, CornerP_good);

    map_P_Vd init_P_IGIGn;
    for (const auto &p : InnerCornerP)
      init_P_IGIGn[p] = {0., 0.};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::vector<map_P_Vd *> P_IGIGn(0);
    for (auto i = 0; i < InnerCornerP.size(); i++) {
      this->P_P_IGIGn[InnerCornerP[i]] = init_P_IGIGn;
      // P_IGIGn.emplace_back(&P_P_IGIGn[InnerCornerP[i]]);
    }

    auto start = std::chrono::high_resolution_clock::now();

#define optimise2

    ///////////////////////////////////////////////////
#ifdef optimise0
    Print("optimise0", Green);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (auto i = 0; i < InnerCornerP.size(); i++) {
      auto &tmp = this->P_P_IGIGn[InnerCornerP[i]];
      for (const auto &f : Faces)
        for (const auto &[p, igign] : BEM::calc_P_IGIGn(f, InnerCornerP[i]))
          tmp[p] += igign;
    }
#endif
    ///////////////////////////////////////////////////
#ifdef optimise1
    Print("optimise1", Green);
    std::vector<map_P_Vd *> tmp(InnerCornerP.size());
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (auto i = 0; i < InnerCornerP.size(); i++)
      tmp[i] = &(this->P_P_IGIGn[InnerCornerP[i]]);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (auto i = 0; i < InnerCornerP.size(); i++)
      for (auto j = 0; j < Faces.size(); j++)
        for (const auto &[p, igign] : BEM::calc_P_IGIGn(Faces[j], InnerCornerP[i]))
          (*tmp[i])[p] += igign;
#endif

    ///////////////////////////////////////////////////
#ifdef optimise2
    Print("optimise2 for deletion of too close ConerP", Green);
    std::vector<map_P_Vd *> tmp(InnerCornerP.size());
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (auto i = 0; i < InnerCornerP.size(); i++)
      tmp[i] = &(this->P_P_IGIGn[InnerCornerP[i]]);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (auto i = 0; i < InnerCornerP.size(); i++)
      for (auto j = 0; j < Faces.size(); j++)
        for (const auto &[p, igign] : BEM::calc_P_IGIGn(Faces[j], InnerCornerP[i]))
          if ((*tmp[i]).find(p) != (*tmp[i]).end())
            (*tmp[i])[p] += igign;
#endif
    ///////////////////////////////////////////////////
#ifdef optimise3
    Print("optimise2 for deletion of too close ConerP", Green);
    std::vector<map_P_Vd *> tmp(InnerCornerP.size());
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (auto i = 0; i < InnerCornerP.size(); i++)
      tmp[i] = &(this->P_P_IGIGn[InnerCornerP[i]]);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (auto i = 0; i < InnerCornerP.size(); i++)
      for (auto j = 0; j < Faces.size(); j++)
        for (const auto &[p, igign] : BEM::calc_P_IGIGn_no_sing_make_sure(Faces[j]->getPoints(), InnerCornerP[i]))
          if ((*tmp[i]).find(p) != (*tmp[i]).end())
            (*tmp[i])[p] += igign;
#endif
    ///////////////////////////////////////////////////

    // for (const auto &p : InnerP)
    //   this->P_P_IGIGn[p][p][1] += 2. * M_PI;
    // for (const auto &p : CornerP)
    //   this->P_P_IGIGn[p][p][1] += M_PI;

    Print("use SolidAngleDN InnerP", Green);
    Print(InnerP.size(), Red);
    for (const auto &p : InnerP) {
      this->P_P_IGIGn[p][p][1] += SolidAngleDN(p);
    }
    Print("use SolidAngleDN CornerP", Green);
    Print(CornerP_good.size(), Red);
    for (const auto &p : CornerP_good) {
      this->P_P_IGIGn[p][p][1] += SolidAngleDN(p);
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << Green << "Elapsed time: " << Red << elapsed.count() << reset << " s\n";

    //=============================================
    //================= exEqs =====================
    //=============================================
    // #define svd1
    map_P_P_Vd exEqs; // 近傍のInnerの数点の平均となるという条件
#if defined(svd1)
    for (const auto &p : CornerP) {
      // 距離に応じて逆数で重みをつける

      // auto neighbors = Intersection(searchPoints(2, p, takeNetworks(InnerCornerP)), InnerCornerP);
      V_netPp neighbors;

      for (auto i = 2; i < 4; i++) {
        neighbors.clear();
        for (const auto &q : searchPoints({i, i}, p, InnerCornerP, takeNetworks(InnerCornerP)))
          if (q->isD() && q != p)
            neighbors.emplace_back(q);
        if (neighbors.size() > 2)
          break;
      }

      // auto neighbors = Intersection(p->getNeighbors(), InnerCornerP);
      V_d r_1({});
      for (auto i = 0; i < neighbors.size(); i++)
        r_1.push_back(1. / Norm(neighbors[i]->getX() - p->getX())); //*影響力*
      double c = 1. / Sum(r_1);

      exEqs[p] = init_P_IGIGn; // 列ベクトルでの初期化
      exEqs[p][p] = {0., 1. /* <- coef of phi*/};

      for (auto i = 0; i < neighbors.size(); i++)
        exEqs[p][neighbors[i]] = {0. /*この方程式にphinは関係ないので0*/, -c * r_1[i]};
    }
    for (auto &[p, p_igign] : exEqs)   // each row
      for (auto &[q, igign] : p_igign) // each column
        if (q->isN())                  /*Neumann*/
          igign = -Reverse(igign);
#endif
    //=============================================
    //==== adjust depending on known or unknown ===
    //=============================================
    // SWAP LEFT AND RIGHT HANDSIDE
    for (auto &[p, p_igign] : this->P_P_IGIGn) // each row
      for (auto &[q, igign] : p_igign)         // each column
        if (q->isN())                          /*Neumann*/
          igign = -Reverse(igign);
    //=============================================
    //================= solve =====================
    //=============================================
    V_d knowns = {};
    Print("(1-2) 境界条件に応じて，matOfUnknowns，matOfKnownsを作成", Green);
    // #ifdef svd
    VV_d igign;
    for (const auto &[org, p_igign] : this->P_P_IGIGn) {
      igign = Transpose(TakeSecond(p_igign));
      this->matOfUnknowns.emplace_back(igign[0]);
      this->matOfKnowns.emplace_back(igign[1]);
      knowns.emplace_back(P_phiphin_BIE[org][org->isD() ? 0 /*phi*/ : 1 /*phin*/]);
    }
    for (const auto &[org, p_igign] : exEqs) {
      igign = Transpose(TakeSecond(p_igign));
      this->matOfUnknowns.emplace_back(igign[0]);
      this->matOfKnowns.emplace_back(igign[1]);
    }

    this->pOfknowns = TakeFirst(this->P_P_IGIGn);

    Print("(2) 未知変数の計算", Green);
    if (this->matOfKnowns.size() == knowns.size())
      solver_name = "lu";
    else
      solver_name = "svd";

    if (solver_name.compare((std::string) "lu") == 0)
      this->lu = new ludcmp(matOfUnknowns /*未知の行列係数（左辺）*/);
    else if (solver_name.compare((std::string) "svd") == 0)
      this->svd = new SVD(matOfUnknowns /*未知の行列係数（左辺）*/);

    V_d phiORphin(knowns.size());

    if (solver_name.compare((std::string) "lu") == 0)
      this->lu->solve(Dot(matOfKnowns, knowns) /*既知のベクトル（右辺）*/, phiORphin /*解*/);
    else if (solver_name.compare((std::string) "svd") == 0)
      this->svd->solve(Dot(matOfKnowns, knowns) /*既知のベクトル（右辺）*/, phiORphin /*解*/);
    else
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");

    Print("use " + solver_name, Red);
    //
    int i = 0;
    for (const auto &p : this->pOfknowns)
      this->P_phiphin_BIE[p][p->isD() ? 1 /*phin*/ : 0 /*phi*/] = phiORphin[i++] /*解*/;
    // BIEの計算終了

    //==============================================================
    //=========== P_phiphin_InnerOuterP (givenを含む) ===============
    //==============================================================
    // V_Netp all_networks = takeNetworks(TakeFirst(P_phiphin_InnerOuterCornerP));

    Print("(ii) OuterPのphiとphinは計算されていないので補間する", Magenta);

    map_P_Vd reference4OuterGivenP = KeyTake(this->P_phiphin_BIE, InnerP); // BIEを解いて得られた結果．これを補間に使用する．
    // map_P_Vd reference4OuterGivenP = this->P_phiphin_BIE; //BIEを解いて得られた結果．これを補間に使用する．

    Print("必ずgiven_P_phiphin_InnerOuterPは含むように", Magenta);
    auto all_OuterGivenP = DeleteDuplicates(Join(OuterP, TakeFirst(given_P_phiphin_InnerOuterP)));

    Print("外部の点のphi,phi_nを補完する", Magenta);
    map_P_Vd P_phiphin_OuterGivenP;

#define outer_point_interpolation
#if defined(outer_point_interpolation)
    ////////////////////////////////////////////////////////// 修正前

    for (const auto &p : all_OuterGivenP)
      P_phiphin_OuterGivenP[p] = {0., 0.};

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < all_OuterGivenP.size(); i++) {
      auto q = all_OuterGivenP[i];
      if (!MemberQ(InnerP, q)) {
        P_phiphin_OuterGivenP[q] = {BEM::IDW_phi(q, reference4OuterGivenP, netsamp.net4phi(q), 3.), BEM::IDW_phin(q, reference4OuterGivenP, netsamp.net4phin(q), 3.)};
      }
    }

    // for (const auto &p : all_OuterGivenP)
    //     if (!MemberQ(InnerP, p))
    //     {
    //         //this->P_phiphin_InnerOuterP[p] = BEM::InterpolationRBF_phiphin(p, fromBIE, {p->getNetwork()});
    //         // P_phiphin_OuterGivenP[p] = BEM::IDW_phiphin(p, fromBIE, getReferenceNetwork(p, all_networks), 6);
    //         auto phi = BEM::IDW_phi(p, reference4OuterGivenP, netsamp.net4phi(p), 3.);
    //         auto phin = BEM::IDW_phin(p, reference4OuterGivenP, netsamp.net4phin(p), 3.); //phinは自身の所属する面だけで補間
    //         P_phiphin_OuterGivenP[p] = {phi, phin};
    //         // P_phiphin_OuterGivenP[p] = BEM::IDW_phiphin(p, reference4OuterGivenP, netsamp.net4phi(p), 3.);
    //         //P_phiphin_OuterGivenP[p] = BEM::InterpolationRBF_phiphin(p, fromBIE, getReferenceNetwork(p, all_networks));
    //     }
#elif defined(outer_point_interpolation_mod)
    //////////////////////////////////////////////////////////////// 修正後
    // この方法だと，角点に近い点が，角に徐々に吸い寄せられる

    // 外部の中で内部の近傍だけをRBFで補間
    // phinも補間するなら，面の法線方向が同じ面を使って補間する必要がある．

    for (const auto &p : all_OuterGivenP)
      P_phiphin_OuterGivenP[p] = {0., 0.};

    V_netPp finishedP({});

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < all_OuterGivenP.size(); i++) {
      auto p = all_OuterGivenP[i];
      if (!MemberQ(InnerP, p)) {
        // 自身の所属するネットワークの格子点は，十分に自分の周りに存在するかどうかチェック
        takeToDepth s({6, 10}, p, reference4OuterGivenP, {p->getNetwork()}, 40);
        // double m = MeanDistance(s.X, p->getX()) / MeanDistance(s.X);
        // std::cout << Red <<"平均的距離 " <<  m << reset<< std::endl;
        // もし十分に存在する場合
        if ((MeanDistance(s.X, p->getX()) / MeanDistance(s.X)) < 1.) {
          auto phi = BEM::InterpolationRBF_phi(p, reference4OuterGivenP, netsamp.net4phi(p), 40);
          auto phin = BEM::InterpolationRBF_phin(p, reference4OuterGivenP, netsamp.net4phin(p), 40); // phinは自身の所属する面だけで補間
          P_phiphin_OuterGivenP[p] = {phi, phin};
#ifdef _OPENMP
#pragma omp critical
#endif
          finishedP.emplace_back(p);
        }
      }
    }
    // RBFで計算した結果をサンプルに追加し，IDWの補間に利用
    reference4OuterGivenP.insert(P_phiphin_OuterGivenP.begin(), P_phiphin_OuterGivenP.end());

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < all_OuterGivenP.size(); i++) {
      auto p = all_OuterGivenP[i];
      if (!MemberQ(InnerP, p) && !MemberQ(finishedP, p)) {
        auto phi = BEM::IDW_phi(p, reference4OuterGivenP, netsamp.net4phi(p), 4.);
        auto phin = BEM::IDW_phin(p, reference4OuterGivenP, netsamp.net4phin(p), 4.); // phinは自身の所属する面だけで補間
        P_phiphin_OuterGivenP[p] = {phi, phin};
      }
    }
#elif defined(outer_point_interpolation_mod_mod)
    //////////////////////////////////////////////////////////////// 修正後
    Print("外部の中で内部の近傍だけをRBFで補間", Magenta);
    Print("phinも補間するなら，面の法線方向が同じ面を使って補間する必要がある．", Magenta);
    V_netPp finishedP({});
    for (const auto &p : all_OuterGivenP)
      if (!MemberQ(InnerP, p)) {
        // 自身の所属するネットワークの格子点は，十分に自分の周りに存在するかどうかチェック
        takeToDepth s({6, 10}, p, reference4OuterGivenP, {p->getNetwork()}, 40);
        double m = MeanDistance(s.X, p->getX()) / MeanDistance(s.X);
        // Print(m,Red);
        // もし十分に存在する場合
        if (m < 1.01) {
          finishedP.emplace_back(p);
          auto phi = BEM::IDW_phi(p, reference4OuterGivenP, netsamp.net4phi(p), 3.);
          auto phin = BEM::IDW_phin(p, reference4OuterGivenP, netsamp.net4phin(p), 3.); // phinは自身の所属する面だけで補間
          P_phiphin_OuterGivenP[p] = {phi, phin};
        }
      }
    Print("RBFで計算した結果をサンプルに追加し，IDWの補間に利用", Magenta);
    reference4OuterGivenP.insert(P_phiphin_OuterGivenP.begin(), P_phiphin_OuterGivenP.end());
    for (const auto &p : all_OuterGivenP)
      if (!MemberQ(InnerP, p) && !MemberQ(finishedP, p)) {
        auto phi = BEM::IDW_phi(p, reference4OuterGivenP, netsamp.net4phi(p), 4.);
        auto phin = BEM::IDW_phin(p, reference4OuterGivenP, netsamp.net4phin(p), 4.); // phinは自身の所属する面だけで補間
        P_phiphin_OuterGivenP[p] = {phi, phin};
        // P_phiphin_OuterGivenP[p] = BEM::IDW_phiphin(p, reference4OuterGivenP, {p->getNetwork()}, 10);
      }
#endif
    /////////////////////////////////////////////////////////////////

    this->P_phiphin_InnerOuterP = KeyTake(this->P_phiphin_BIE, InnerP);
    auto tmpmap = P_phiphin_OuterGivenP;
    this->P_phiphin_InnerOuterP.Collapse(tmpmap);
    //==============================================================
    //========= P_phiphin_Inne∞rOuterCornerP (givenを含む) ==========
    //==============================================================
    this->P_phiphin_InnerOuterCornerP = this->P_phiphin_InnerOuterP;
    tmpmap = this->P_phiphin_BIE;
    this->P_phiphin_InnerOuterCornerP.Collapse(tmpmap);
    //==============================================================
    //========================= nablaPhi ===========================
    //==============================================================
    Print("(iii) nablaPhiを計算（外部の流速も計算する）", Magenta);

    // uo_map_P_Vd uomap_P_phiphin_InnerOuterCornerP = map2unordered_map(P_phiphin_InnerOuterCornerP);
    // uo_map_P_Vd uomap_P_phiphin_InnerOuterP = map2unordered_map(P_phiphin_InnerOuterP);

#define _do_not_use_cornerP
#if defined(_use_RBF_for_all_mod)
    /*内部はBIEを使う*/
    // 全てRBF補間の微分で流速を評価
    for (const auto &p : TakeFirst(P_phiphin_InnerOuterCornerP)) {
      if (MemberQ(all_OuterGivenP, p))
        this->P_nablaPhi[p] = BEM::InterpolationRBFvG(p, P_phiphin_InnerOuterCornerP, getReferenceNetwork(p, all_networks));
      else
        this->P_nablaPhi[p] = BEM::InterpolationRBFvG(p, P_phiphin_BIE, getReferenceNetwork(p, all_networks));
    }
#elif defined(_do_not_use_cornerP)
    /*内部はBIEを使う*/
    // 全てRBF補間の微分で流速を評価

    auto ps = TakeFirst(P_phiphin_InnerOuterCornerP);
    for (const auto &p : ps)
      this->P_nablaPhi[p] = {0., 0., 0.};

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (auto i = 0; i < ps.size(); i++) {
      if (MemberQ(CornerP, ps[i]))
        this->P_nablaPhi[ps[i]] = BEM::InterpolationRBFvG(ps[i], P_phiphin_InnerOuterCornerP, netsamp.net4nabla(ps[i]));
      else
        this->P_nablaPhi[ps[i]] = BEM::InterpolationRBFvG(ps[i], P_phiphin_InnerOuterP, netsamp.net4nabla(ps[i]));
    }
    // for (const auto &p : TakeFirst(P_phiphin_InnerOuterCornerP))
    // {
    //     if (MemberQ(CornerP, p))
    //         this->P_nablaPhi[p] = BEM::InterpolationRBFvG(p, uomap_P_phiphin_InnerOuterCornerP, netsamp.net4nabla(p));
    //     else
    //         this->P_nablaPhi[p] = BEM::InterpolationRBFvG(p, uomap_P_phiphin_InnerOuterP, netsamp.net4nabla(p));
    // }
#elif defined(_use_RBF_for_all_)
    // 全てRBF補間の微分で流速を評価
    for (const auto &p : TakeFirst(P_phiphin_InnerOuterCornerP)) {
      // this->P_nablaPhi[p] = BEM::InterpolationRBFvG(p, P_phiphin_InnerOuterCornerP, (p->isC() ? takeN(all_networks) : V_Netp{p->getNetwork()}));

      this->P_nablaPhi[p] = BEM::InterpolationRBFvG(p, P_phiphin_InnerOuterCornerP, getReferenceNetwork(p, all_networks));

      // this->P_nablaPhi[p] = BEM::meanvG(p,P_phiphin_InnerOuterCornerP,(p->isC() ? takeN(all_networks) : V_Netp{p->getNetwork()}),cpg.well_Faces);
      //  this->P_nablaPhi[p] = BEM::global_velocity(p,
      //                                             P_phiphin_InnerOuterCornerP,
      //                                             (p->isC() ? takeN(all_networks) : V_Netp{p->getNetwork()}),
      //                                             cpg.well_Faces);
    }
#elif defined(_use_IDW_for_all_)
    for (const auto &p : TakeFirst(P_phiphin_InnerOuterCornerP)) {
      if (MemberQ(all_OuterGivenP, p))
        this->P_nablaPhi[p] = BEM::InterpolationIDWvG(p, P_phiphin_InnerOuterCornerP, getReferenceNetwork(p, all_networks), 2.);
      else
        this->P_nablaPhi[p] = BEM::InterpolationRBFvG(p, P_phiphin_InnerOuterCornerP, getReferenceNetwork(p, all_networks));
    }
#elif defined(_use_RBF_for_all_step_by_step_)
    V_netPp InnerP_({}), OuterP_({}), CornerP_({}), OtherP_({});
    for (const auto &p : TakeFirst(P_phiphin_InnerOuterCornerP)) {
      if (MemberQ(InnerP, p))
        InnerP_.emplace_back(p);
      if (MemberQ(OuterP, p))
        OuterP_.emplace_back(p);
      if (MemberQ(CornerP, p))
        CornerP_.emplace_back(p);
      else
        OtherP_.emplace_back(p);
    }

    for (const auto &p : InnerP_)
      this->P_nablaPhi[p] = BEM::InterpolationRBFvG(p, P_phiphin_InnerOuterCornerP, getReferenceNetwork(p, all_networks));

    auto ref = this->P_nablaPhi;
    for (const auto &p : OuterP_) {
      takeToDepth tmp({3, 10}, p, ref, getReferenceNetwork(p, all_networks));
      InterpolationIDW interp(tmp.X, tmp.VV, 5.);
      this->P_nablaPhi[p] = interp.calcVV(p->getX());
    }
    for (const auto &p : CornerP_)
      this->P_nablaPhi[p] = BEM::InterpolationRBFvG(p, P_phiphin_InnerOuterCornerP, getReferenceNetwork(p, all_networks));

    for (const auto &p : OtherP_)
      this->P_nablaPhi[p] = BEM::InterpolationRBFvG(p, P_phiphin_InnerOuterCornerP, getReferenceNetwork(p, all_networks));

#endif

    //==============================================================
    //=========== BIEで計算しなかった CornerP_bad の補間 ===============
    //==============================================================
    //         for (const auto &q : CornerP_bad)
    //         {
    //             this->P_phiphin_BIE[q] = {BEM::InterpolationRBF_phi(q, this->P_phiphin_BIE, netsamp.net4phi(q)),
    //                                       BEM::InterpolationRBF_phin(q, this->P_phiphin_BIE, netsamp.net4phin(q))};
    //         }

    // #ifdef _OPENMP
    // #pragma omp parallel for
    // #endif
    //         for (auto i = 0; i < CornerP_bad.size(); i++)
    //         {
    //             this->P_nablaPhi[CornerP_bad[i]] = BEM::InterpolationRBFvG(CornerP_bad[i], P_phiphin_InnerOuterCornerP, netsamp.net4nabla(ps[i]));
    //         }

    //==============================================================
    //========================= P_DphiDt ===========================
    //==============================================================
    Print("(5) DphiDtをベルヌーイの式から計算．流速の計算には放射関数基底補間の微分を用いる．（外部の流速も計算する）", Green);
    for (const auto &[p, u] : P_nablaPhi)
      this->P_DphiDt[p] = Dot(u, u) / 2. - p->getX()[2]; // これは圧力pが含まれていないので，Neumann境界条件では正しくない
  };

  //---------------------------------------------------------------

  /*BIEを使って，D->N,N->Dを実行する．既知の変数を引数として入力する．もし，与えられなかった変数があれば，自動で0が代入されBIEを解く*/
  /*NeumannとDirichletの混同した結果が帰ってくるので注意*/
  map_P_Vd knowns_DN;
  map_P_d solveBIE(map_P_d knowns_D_IN, map_P_d knowns_N_IN = {}) {
    V_d knowns(0);
    // とりあえずknowns_INを入れてもいい．わからないものは，0を代入する．
    // NeumannならphinなどDirichletならphiなどがknownsに代入される．
    Print("NeumannならphinなどDirichletならphiなどがknownsに代入される．", red);
    for (const auto &p : this->pOfknowns) {
      map_P_d::iterator it;
      if ((it = knowns_D_IN.find(p)) != knowns_D_IN.end())
        knowns.emplace_back(it->second);
      else if ((it = knowns_N_IN.find(p)) != knowns_N_IN.end())
        knowns.emplace_back(it->second);
      else
        knowns.emplace_back(0);
    }
    //
    V_d solution(knowns.size());

    if (!(solver_name.compare((std::string) "lu") == 0))
      this->lu->solve(Dot(this->matOfKnowns, knowns) /*既知のベクトル（右辺）*/, solution /*解*/);
    else if (solver_name.compare((std::string) "svd") == 0)
      this->svd->solve(Dot(this->matOfKnowns, knowns) /*既知のベクトル（右辺）*/, solution /*解*/);

    //
    int i = 0;
    map_P_d ret;
    for (const auto &p : this->pOfknowns)
      ret[p] = solution[i++];

    Print("solutionはDirichletとNeumannが混合しているので扱いにくい．そこでN,Dで整理したベクトルを作成しメンバとして保存しておく", red);
    knowns_DN.clear();
    for (const auto &p : this->pOfknowns) {
      map_P_d::iterator it;
      if (p->isD()) {
        if ((it = knowns_D_IN.find(p)) != knowns_D_IN.end())
          knowns_DN[p] = V_d{it->second, ret[p]};
        else
          knowns_DN[p] = V_d{0., ret[p]};
      } else if (p->isN()) {
        if ((it = knowns_N_IN.find(p)) != knowns_N_IN.end())
          knowns_DN[p] = V_d{ret[p], it->second};
        else
          knowns_DN[p] = V_d{ret[p], 0.};
      } else
        knowns_DN[p] = V_d{0., 0.};
    }

    return ret;
  };
};

#endif