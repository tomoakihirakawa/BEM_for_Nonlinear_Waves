#ifndef CompGrid_H
#define CompGrid_H
#pragma once

/*CompGrid_detail
複数の`selector`の格子をまとめあげ，計算格子を作る．
メンバ変数として，少し変わった干渉ネットワーク`xnet`を作成する．
`xnet`は，自身のメンバ変数として面を作り内部に保持するが，
面のnetworkは`xnet`を指すポインタではなく，`xnet`が対象とした2つのネットワークのどちらかを指す（面が乗っているどちらかのネットワークを指している）．
しかし，面を構成する辺のネットワークは，`xnet`を指すポインタである．このようにする理由は，面の境界条件をわかりやすくするためである．
`CompGrid`オブジェクトのデストラクタは，`xnet`を`delte`する．
CompGrid_detail*/
/*CompGrid_code*/
class CompGrid {
private:
  //-----------------------------------------
  bool isWellCondition(const netFp &f) const {
    for (const auto &ang : f->getAngles())
      if (ang < 1E-3)
        return false;
    if (f->getArea() < 1E-6)
      return false;
    return true;
  };
  //-----------------------------------------
  // bool testGetNormal(const netPp &p) const
  // {
  //     ////////////// get normal test //////////
  //     VV_d normals({});
  //     V_netFp fs;
  //     V_netPp ps;
  //     for (auto d = 1; d <= 2; d++)
  //     {
  //         //------------------
  //         if (d == 1)
  //             fs = p->getBoundaryFaces(); //これだけならgetNormal()と同じ
  //         else
  //         {
  //             depth_searcher<networkPoint> S(d);
  //             S.set(p);
  //             S.search();
  //             ps = S.getObjects();
  //             fs.clear();
  //             for (const auto &p : ps)
  //                 for (const auto &f : p->getBoundaryFaces())
  //                     fs.emplace_back(f);
  //             fs = DeleteDuplicates(fs);
  //         }
  //         //------------------
  //         if (p->isD())
  //         {
  //             for (const auto &f : takeD(fs))
  //                 normals.emplace_back(f->getNormal());
  //         }
  //         else
  //         {
  //             for (const auto &f : takeN(fs))
  //                 normals.emplace_back(f->getNormal());
  //         }
  //         //------------------
  //         if (!normals.empty())
  //             return true;
  //     }
  //     ///////////////////
  //     return false;
  // };
  //-----------------------------------------
  // bool isWellCondition(const netPp &p) const
  // {
  //     for (const auto &f : p->getBoundaryFaces())
  //         if (!isWellCondition(f))
  //             return false;
  //     return testGetNormal(p);
  // };

public:
  V_netPp mandatoryP;
  V_netPp InnerP, CornerP, OuterP, PointsForBIE, InnerOuterP, InnerCornerP;
  V_netFp all_Faces, well_Faces, ill_Faces;
  VV_netPp routeP;
  Network *xnet;
  // std::vector<selector *> s0;
  // std::vector<selector *> s1;
  // BEM::DN nets;

  ~CompGrid() { delete (this->xnet); };

  CompGrid() : mandatoryP({}), InnerP({}), CornerP({}), OuterP({}), PointsForBIE({}), InnerOuterP({}), InnerCornerP({}), all_Faces({}), well_Faces({}), ill_Faces({}), routeP({}), xnet(nullptr) {};

  CompGrid(const V_Netp &base_net, const VV_d &cardinal, const V_netPp mandatoryP_IN = {}) : mandatoryP({}), InnerP({}), CornerP({}), OuterP({}), PointsForBIE({}), InnerOuterP({}), InnerCornerP({}), all_Faces({}), well_Faces({}), ill_Faces({}), routeP({}), xnet(nullptr) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    this->xnet = new Network(base_net, "Xnetwork");

    mk_vtu("./vtu/cpg_" + xnet->getName() + "_Lines_init.vtu", this->xnet->getLines());

    calcCompGrid(base_net, cardinal);

    display(this->xnet);

    mk_vtu("./vtu/cpg_" + xnet->getName() + "_Lines.vtu", this->xnet->getLines());
    mk_vtu("./vtu/cpg_" + xnet->getName() + "_Faces.vtu", this->xnet->getBoundaryFaces());
    mk_vtu("./vtu/cpg_" + xnet->getName() + "_all_Faces.vtu", this->all_Faces);
    mk_vtu("./vtu/cpg_" + xnet->getName() + "_routeP.vtu", this->routeP);
    mk_vtu("./vtu/cpg_" + xnet->getName() + "_Lines_Intxn.vtu", takeIntxn(this->xnet->getLines()));
  };
  ////////////////////////////////
  bool canFlip2(const netLp l, const V_Netp &target_storage_IN) {
    if (l->isIntxn())
      return false;

    auto fs = l->getBoundaryFaces();
    if (fs.size() != 2)
      return false;

    auto anet = fs[0]->getNetwork();
    for (const auto &f : fs) {
      if (anet != f->getNetwork())
        return false;
      if (!(MemberQ(target_storage_IN, f->getStorage())))
        return false;
    }
    return true;
  };
  bool canFlip(const netLp l, const V_Netp &target_storage_IN) {
    auto fs = l->getBoundaryFaces();
    if (isFlat(l, 1E-2) && fs.size() == 2 && !l->isIntxn() && (fs[0]->getStorage() == fs[1]->getStorage()) && (fs[0]->getNetwork() == fs[1]->getNetwork()))
      if (MemberQ(target_storage_IN, fs[0]->getStorage()))
        if (MemberQ(target_storage_IN, fs[1]->getStorage()))
          return true;
    return false;
  };

  class isVeryBadCase {
  public:
    V_netFp fs;
    V_netLp ls;
    isVeryBadCase() : fs({}), ls({}) {};
    V_netLp operator()(const netLp l, const double minarea) {
      fs.clear();
      ls.clear();
      V_netLp next_l = {l};
      do {
        auto ll = next_l;
        next_l.clear();
        for (const auto &l : ll)
          for (const auto &f : l->getBoundaryFaces())
            if (!l->isIntxn() && !MemberQ(fs, f) && f->getArea() < minarea) {
              fs.emplace_back(f);
              ls.emplace_back(l);
              next_l.emplace_back(l);
            }
      } while (!next_l.empty());
      return ls;
    };
  };

  /////////////////////////////////
  void remesh_flipIfIllegal(Network *net) {
    try {
      std::cout << __PRETTY_FUNCTION__ << std::endl;
      auto PS = net->getPoints();
      LaplacianSmoothing_IfOnStraight_IntXLine(PS /*内部でシャッフルする*/); // ちゃんとsetXされているかチェック
      bool found = false;
      int c = 0, mergecount = 0;
      do {

        LaplacianSmoothing_IfOnStraight_IntXLine(net->getPoints());

        found = false;
        auto lines = net->getLines();
        std::shuffle(std::begin(lines), std::end(lines), std::default_random_engine());
        for (const auto &l : lines) {
          // if (canFlip(l, {net}))
          //     found = l->flipIfIllegal();

          if (canFlip(l, {net})) {
            found = l->flipIfIllegal();
            auto fs = l->getBoundaryFaces();
            for (auto i = 0; i < 2; i++)
              for (auto &f : fs)
                LaplacianSmoothing_IfOnStraight_IntXLine(f->getPoints() /*内部でシャッフルする*/); // ちゃんとsetXされているかチェック
          }
        }

        LaplacianSmoothing_IfOnStraight_IntXLine(net->getPoints());

        lines = net->getLines();
        std::shuffle(std::begin(lines), std::end(lines), std::default_random_engine());

        double r;
        for (const auto &l : lines)
          if (canFlip2(l, {net})) {
            auto fs = l->getBoundaryFaces();
            r = fs[0]->getArea() / fs[1]->getArea();
            if (r < 1E-5 || r > 1E+5) {
              std::cout << "r = " << r;
              found = l->Flip();
              for (auto i = 0; i < 2; i++)
                for (auto &f : fs)
                  LaplacianSmoothing_IfOnStraight_IntXLine(f->getPoints() /*内部でシャッフルする*/); // ちゃんとsetXされているかチェック
              std::cout << Red << " -> " << fs[0]->getArea() / fs[1]->getArea() << reset << std::endl;
            }
          }

        // for (const auto &l : lines)
        // {
        //     if (canFlip2(l, {net}))
        //     {
        //         auto fs = l->getBoundaryFaces();
        //         r = Min(V_d{fs[0]->getArea() / fs[1]->getArea(),fs[1]->getArea() / fs[0]->getArea()});
        //         if (r < 1E-3)
        //         {
        //             auto current_r = r;
        //             std::cout << "r = " << r << std::endl;
        //             if(l->Flip())
        //             {
        //                 auto _fs = l->getBoundaryFaces();
        //                 if(current_r > Min(V_d{_fs[0]->getArea() / _fs[1]->getArea(),_fs[1]->getArea() / _fs[0]->getArea()}))
        //                     l->Flip();
        //             }
        //         }
        //     }
        // }

        // for (const auto &L : lines)
        // {
        //     isVeryBadCase vbc;
        //     auto ls = vbc(L, 1E-5);
        //     Print(ls);
        //     for (const auto &l : ls)
        //     {
        //         if (canFlip2(l, {net}))
        //         {
        //             for (const auto &f : l->getBoundaryFaces())
        //             {
        //                 auto ps = f->getPoints(l);
        //                 if (!ps[2]->isXPoint())
        //                 {
        //                     found = l->flipIfIllegal();
        //                 }
        //             }
        //         }
        //     }
        // }

        // for (auto i = 0; i < 100; i++)
        // {
        //     for (const auto &l : net->getLines())
        //     {
        //         if (l->isIntxn() && l->length() < 1E-2)
        //             if (l->getPoints()[0]->getNetwork() == net && l->getPoints()[1]->getNetwork() == net)
        //             {
        //                 std::cout << "l->length() = " << l->length() << std::endl;
        //                 if (l->Collapse())
        //                 {
        //                     mergecount++;
        //                     for (auto &f : l->getBoundaryFaces())
        //                         LaplacianSmoothing_IfOnStraight_IntXLine(f->getPoints() /*内部でシャッフルする*/); //ちゃんとsetXされているかチェック
        //                     break;
        //                 }
        //                 // std::cin.ignore();
        //             }
        //     }
        // }

        for (auto i = 0; i < 100; i++) {
          for (const auto &p : net->getPoints()) {
            V_netPp ps({});
            V_netLp intxlines({});
            netPp q;
            for (const auto &l : p->getLines())
              if (l->isIntxn())
                if ((q = (*l)(p)) && q->isXPoint() && p->isXPoint()) {
                  network::add(ps, q);
                  network::add(intxlines, l);
                }

            if (ps.size() == 2) {
              // Print("平滑化可能");
              auto v0 = p->getX() - ps[0]->getX();
              auto v1 = ps[1]->getX() - p->getX();
              auto angle = MyVectorAngle(v0, v1);
              if (isFinite(angle) && std::abs(angle) < 1E-2) {
                if (intxlines[0]->length() < 1E-2) {
                  intxlines[0]->Collapse();
                  break;
                }
              }
            }
          }
        }

        // ここを変更したがうまくいかない．．
        std::cout << "c = " << c++ << std::endl;
      } while (c <= 5 /*最低回数*/ || (found && c < 50));
      std::cout << "mergecount = " << mergecount << std::endl;
    } catch (const error_message &e) {
      throw;
    }
  };
  ////////////////////////////////////////////////
  std::map<Network *, V_netFp> takePossibleXFaces(const V_Netp &base_net, const VV_d &cardinals) {
    Print("三角分割する可能性がある面を抜き出す", lRed);
    std::map<Network *, V_netFp> ret; //{{&water,water_XFaces}, {&obj,obj_XFaces}};
    try {
      for (const auto &cardinal : cardinals)
        for (const auto &net : base_net)
          if (net->isB()) {
            Print("0番は，交線と接続している面を抜き出す", Red);
            selector s({this->xnet}); //<==================================================
            s.takePossibleXFaces(net, cardinal);
            network::add(ret[net], s.faces);
#if defined(simulation)
            mk_vtu("./vtu/" + net->getName() + "_kanousei" + std::to_string(time_step) + ".vtu", extractPoints(s.faces));
#endif
          }
      return ret;
    } catch (const error_message &e) {
      throw;
    }
  };
  //------------------------------------------------------
  void makeXFaces(const V_Netp &base_net, const VV_d &cardinals) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    try {
      ///////// しょうがなく
      // std::map<Network *, V_netFp> map;
      // for (const auto &n : base_net)
      // {
      //     map[n] = n->getBoundaryFaces();
      //     mk_vtu("./vtu/makeFaces_" + n->getName() + "->Faces" + std::to_string(time_step) + ".vtu", n->getBoundaryFaces());
      // }
      // std::cin.ignore();
      ////
      // Print(message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "XNetwork.Facesの作成"));
      // for (const auto &[netp, fs] : takePossibleXFaces(base_net, cardinals))

      // for (const auto &[netp, fs] : map)
      for (const auto &[netp, fs] : takePossibleXFaces(base_net, cardinals))
        for (const auto &f : fs)
          for (const auto &polygon : f->getPointsCutFaces()) {
            this->routeP.emplace_back(polygon);
            if (polygon.size() >= 3) {
              auto V_trp_ps = triangulate(polygon, f->getNormal(), 1E-8 /*むしろこの値が大きいことがエラーにつながる*/);
              for (const auto &tri_ps : V_trp_ps)
                new networkFace(netp /*所属*/, this->xnet /*保存先*/, link(tri_ps, this->xnet /*所属保存先*/)); // 変更2021/01/27コンストラクタで自動保存&自動保存先保存. 作られた面と与えられた面だけが繋がっている状態
            } else
              throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "polygon size is < 3");
          }
    } catch (const error_message &e) {
      throw;
    }
    Print("makeFaces done");
  };
  //------------------------------------------------------
  V_netFp takeRelatedFaces(const V_Netp &base_nets, const VV_d &cardinals) {
    V_netFp ret;
    Print("relatedFacesの抽出，", lRed);
    for (auto i = 0; i < base_nets.size(); i++) {
      auto c = cardinals[(i < cardinals.size() ? i : cardinals.size() - 1)];
      if (base_nets[i]->isB()) {
        selector s({this->xnet});
        s.takeFaces(base_nets[i], c);
        for (const auto &f : s.faces)
          ret.emplace_back(f);
      }
    }
    return ret;
  };
  //------------------------------------------------------
  template <typename T> void Delete(std::vector<T *> ps) {
    auto tmpFs = DeleteDuplicates(ps); // important
    for (auto i = 0; i < tmpFs.size(); i++) {
      tmpFs[i]->Delete();
      delete tmpFs[i];
      tmpFs[i] = nullptr;
    }
  };

  void deleteUnrelatedFaces(const V_Netp &base_nets, const VV_d &cardinals) {
    V_netFp relatedFaces = takeRelatedFaces(base_nets, cardinals), unrelatedFaces = {};
    for (const auto &f : DeleteDuplicates(this->xnet->getBoundaryFaces()))
      if (!MemberQ(relatedFaces, f))
        unrelatedFaces.emplace_back(f);

    mk_vtu("./vtu/relatedFaces" + std::to_string(time_step) + ".vtu", relatedFaces);
    mk_vtu("./vtu/unrelatedFaces" + std::to_string(time_step) + ".vtu", unrelatedFaces);

    Print("無関係の面を削除", Red);
    for (auto i = 0; i < unrelatedFaces.size(); i++)
      delete unrelatedFaces[i];

    Print("無関係の点を削除", Red);
    auto tmpPs = DeleteDuplicates(this->xnet->getPoints()); // important
    for (auto i = 0; i < tmpPs.size(); i++)
      if (tmpPs[i]->getBoundaryFaces().empty())
        delete tmpPs[i];
    Print("deleteUnrelatedFaces done");
  };
  //--------------------------------------------------------
  // V_netFp takeIllFaces(const V_netFp &fs)
  // {
  //     V_netFp ret;
  //     for (const auto &f : fs)
  //         if (!isWellCondition(f))
  //             ret.emplace_back(f);
  //     return ret;
  // };
  // V_netFp takeWellFaces(const V_netFp &fs)
  // {
  //     V_netFp ret;
  //     for (const auto &f : fs)
  //         if (isWellCondition(f))
  //             ret.emplace_back(f);
  //     return ret;
  // };
  //-------
  void calcCompGrid(const V_Netp &base_net, const V_d &cardinal) { calcCompGrid(base_net, {cardinal}); };
  void calcCompGrid(const V_Netp &base_net, const VV_d &cardinals) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    /////////////////
    auto start = std::chrono::high_resolution_clock::now();
    //////////////////
    /*三角分割によって干渉ネットワークに面を挿入する前段階として，三角分割する可能性がある面を抜き出し，三角分割する．*/
    makeXFaces(base_net, cardinals);
    deleteUnrelatedFaces(base_net, cardinals);
    remesh_flipIfIllegal(this->xnet);

    Print("InnerP, CornerP, OuterP, Facesの取得", lRed);
    for (const auto &cardinal : cardinals)
      for (const auto &net : base_net)
        if (net->isB()) {
          selector s({this->xnet}); //<==================================================
          // this->s1.emplace_back(s);
          s.takeFaces(net, cardinal);

          this->InnerP.insert(this->InnerP.end(), s.innerP.begin(), s.innerP.end());
          this->CornerP.insert(this->CornerP.end(), s.cornerP.begin(), s.cornerP.end());
          this->OuterP.insert(this->OuterP.end(), s.outerP.begin(), s.outerP.end());

          // この時点でまだ，s->facesに十分に点が含まれていない
          for (const auto &f : s.faces) {
            this->all_Faces.emplace_back(f);
            if (isWellCondition(f))
              this->well_Faces.emplace_back(f);
            else
              this->ill_Faces.emplace_back(f);
          }
          s.show();
        }

    Print("DeleteDuplicates", Red);
    this->well_Faces = DeleteDuplicates(well_Faces);
    this->ill_Faces = DeleteDuplicates(ill_Faces);
    this->all_Faces = DeleteDuplicates(all_Faces);
    this->InnerP = DeleteDuplicates(InnerP);
    this->CornerP = DeleteDuplicates(CornerP);
    this->OuterP = DeleteDuplicates(OuterP);
    //
    Print("k回Outer拡張 重要：!q->isXPoint()", Red);
    for (auto k = 0; k < 1; k++) {
      V_netPp exOuterP;
      for (const auto &p : OuterP)
        for (const auto &q : p->getNeighbors())
          if (!q->isXPoint() && !MemberQ(CornerP, q) && !MemberQ(InnerP, q) && !MemberQ(OuterP, q) && !MemberQ(exOuterP, q))
            exOuterP.emplace_back(q);

      this->OuterP.insert(this->OuterP.end(), exOuterP.begin(), exOuterP.end());
    }
    // auto tmp = this->mandatoryP;
    // network::erase(tmp, this->InnerP);
    auto tmp = TakeExcept(this->mandatoryP, this->InnerP);
    this->OuterP = DeleteDuplicates(Join(OuterP, tmp));
    this->InnerOuterP = DeleteDuplicates(Join(this->InnerP, this->OuterP));
    this->InnerCornerP = DeleteDuplicates(Join(this->InnerP, this->CornerP));

    this->PointsForBIE.clear();
    Print("PointsForBIEの取得", Red);
    this->PointsForBIE.insert(PointsForBIE.end(), InnerP.begin(), InnerP.end());

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << Green << "CompGrid Elapsed time: " << Red << elapsed.count() << reset << " s\n";
    ///////////////////////////
  };
};
/*CompGrid_code*/

#endif