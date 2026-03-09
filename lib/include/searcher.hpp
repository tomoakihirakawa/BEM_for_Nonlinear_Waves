#ifndef searcher_H
#define searcher_H

#include "Network.hpp"

/*searcher_detail
`searcher`は，

- **`condEnterLine`**の条件にしたがって，`networkLine`を進み，
- **`condGetObject`**の条件にしたがって`networkObject`を取得していく．
- **`condKeepSearch`**の条件にしたがって，探査を続けるかどうか決定する

**また，`searcher`は条件に関わらず，一度通った`networkLine`を再び通ることない．**
`networkObject`クラスのdoUKM, remember, forget関数を使って，通過した`searcher`を管理している．
searcher_detail*/
/*searcher_code*/
/*   @-->@-->@   */
/*   |   |       */
/*   V   |       */
/*   @   @-->@   */
#define use_binary_search
template <class T> class searcherCommon {
protected:
  V_netLp enteredLines, enteredLines_, reachedLines, reachedLines_;
  V_netLp penetrateLines;   // 特別
  std::vector<T *> netObjs; /*取得したオブジェクト*/
  std::vector<T *> netObjs_ /*取得しなかったオブジェクト*/;
  std::vector<T *> netObjs__; /*特別な方法で取得したオブジェクト*/

  // T *startObj;
  // bool TorF;
  std::vector<Network *> networks;
  Network *network;
  /*深さ毎に取得したオブジェクト．[0]は探査をスタートするオブジェクトなので常に1つ*/
  std::vector<std::vector<T *>> netObjsAtDepth;
  /* ------------------------------------------------------ */
  std::unordered_set<networkLine *> enteredLinesUO, enteredLines_UO, reachedLinesUO, reachedLines_UO;
  std::unordered_set<T *> netObjsUO;                     /*取得したオブジェクト*/
  std::unordered_set<T *> netObjs_UO;                    /*取得しなかったオブジェクト*/
  std::unordered_set<T *> netObjs__UO;                   /*特別な方法で取得したオブジェクト*/
  std::vector<std::unordered_set<T *>> netObjsAtDepthUO; /*深さ毎に取得したオブジェクト．[0]は探査をスタートするオブジェクトなので常に1つ*/

public:
  ////////// 並列化のための新しい関数 //////////////////
  // bool doIKU(const netLp l) const
  // {
  // 	return (std::binary_search(this->enteredLines.cbegin(), this->enteredLines.cend(), l) ||
  // 			std::binary_search(this->enteredLines_.cbegin(), this->enteredLines_.cend(), l));
  // 	// return (std::find(this->enteredLines.cbegin(), this->enteredLines.cend(), l) != list.cend())
  // 	// MemberQ(this->enteredLines, l) || MemberQ(this->enteredLines_, l));
  // };
  // bool doIKU(const T *p) const
  // {
  // 	return (std::binary_search(this->netObjs.cbegin(), this->netObjs.cend(), p) ||
  // 			std::binary_search(this->netObjs_.cbegin(), this->netObjs_.cend(), p));
  // 	// return (MemberQ(this->netObjs, p) || MemberQ(this->netObjs_, p) || MemberQ(this->netObjs__, p));
  // };

  //////////////////////////////////////////////////
  int depth;
  //-------------------
  const V_netLp &getEnteredLines() const { return enteredLines; };
  const V_netLp &getEnteredLines_() const { return enteredLines_; };
  const V_netLp &getReachedLines() const { return reachedLines; };
  const V_netLp &getReachedLines_() const { return reachedLines_; };
  const V_netLp &getPenetrateLines() const { return penetrateLines; }; // 特別
  bool addNetwork(Network *net) { return network::add(networks, net); };
  std::vector<bool> addNetworks(const std::vector<Network *> &net) {
    std::vector<bool> ret({});
    for (const auto &n : net)
      ret.push_back(network::add(networks, n));
    return ret;
  };
  bool eraseNetwork(Network *net) { return network::erase(networks, net); };
  std::vector<Network *> getNetworks() { return networks; };

  void setNetwork(Network *net) { this->network = net; };
  void setNetworks(const std::vector<Network *> &nets) { this->networks = nets; };
  Network *getNetwork(Network *net) const { return this->network; };
  const std::vector<T *> &getObjects() const { return this->netObjs; };
  const std::vector<T *> &getObjects_() const { return this->netObjs_; };
  const std::vector<T *> &getObjects__() const { return this->netObjs__; };
  const std::vector<T *> &getAllObjects() const { return Flatten(std::vector<std::vector<T *>>{this->netObjs, this->netObjs_, this->netObjs__}); };
  const std::vector<std::vector<T *>> &getObjectsAtDepth() const { return this->netObjsAtDepth; };

  const std::unordered_set<T *> &getObjectsUO() const { return this->netObjsUO; };
  const std::unordered_set<T *> &getObjects_UO() const { return this->netObjs_UO; };
  const std::unordered_set<T *> &getObjects__UO() const { return this->netObjs__UO; };
  // std::unordered_set<T *> getAllObjectsUO() const { return Flatten(std::vector<std::vector<T *>>{this->netObjsUO, this->netObjs_UO, this->netObjs__UO}); };
  const std::vector<std::unordered_set<T *>> &getObjectsAtDepthUO() const { return this->netObjsAtDepthUO; };
  void clear() {
    enteredLines.clear();
    enteredLines_.clear();
    reachedLines.clear();
    reachedLines_.clear();
    penetrateLines.clear();
    netObjs.clear();
    netObjs_.clear();
    netObjs__.clear();
    networks.clear();
    // startObj = nullptr;
    network = nullptr;
    netObjsAtDepth.clear();
    netObjsAtDepthUO.clear();
    /* ------------------------------------------------------ */
    enteredLinesUO.clear();
    enteredLines_UO.clear();
    reachedLinesUO.clear();
    reachedLines_UO.clear();
    netObjsUO.clear();
    netObjs_UO.clear();
    netObjs__UO.clear();
  };
  void set(T *start_) {
    clear();
    netObjs = {start_};
    netObjsUO = {start_};
    // startObj = (start_);
    network::add(networks, start_->getNetwork());
  };
  void set(std::vector<T *> start_) {
    clear();
    std::sort(start_.begin(), start_.end());
    netObjs = start_;
    netObjsUO.insert(start_.begin(), start_.end());
    // startObj = start_;
    for (const auto &o : start_)
      network::add(networks, o->getNetwork());
  };

  searcherCommon() : enteredLines({}), enteredLines_({}), reachedLines({}), reachedLines_({}), penetrateLines({}), netObjs({}), netObjs_({}), netObjs__({}), /*startObj(nullptr), */ networks({}), network(nullptr), netObjsAtDepth({}) {};

  ///*condEnterLine_*/`condEnterLine`は，`searcher`が目の前にある`networkLine`を通るかどうかを決定する．判断材料は現在地の`networkObect`と目の前にある`networkLine`の持つ内容．**一度チェックした`networkLine`は二度と通らないので注意．**/*condEnterLine_*/
  virtual bool condEnterLine(const T *p, const netLp l) { return true; };
  ///*condGetObject_*/`condGetObject`は，`searcher`が現在地の（手元にある）`networkObject`を取得するかどうかを決定する．判断材料は通ってしまい後ろにある`networkLine`と現在地の`networkObject`の持つ内容．/*condGetObject_*/
  virtual bool condGetObject(const netLp l, const T *P) { return true; };
  ///*condKeepSearch_*/`condKeepSearch`は，`searcher`が現在地の（手元にある）`networkObject`から，再び探査を開始するかどうかを判断する．判断材料は通ってしまい後ろにある`networkLine`と現在地の`networkObject`の持つ内容．/*condKeepSearch_*/
  virtual bool condKeepSearch(const netLp l, const T *P) { return true; };
  virtual void initSearch() {};
};

/////////////////////////////////////////////////////////////////////////////
template <typename T> class searcher;

template <> class searcher<networkPoint> : public searcherCommon<networkPoint> {
public:
  searcher() : searcherCommon<networkPoint>() {};

  /*search_detail
点と点は，直接繋がっておらず，間には`networkLine`が存在する．これによって，`search()`は，`condEnterLine`を使って，通過する線の条件も変更することができる．このことは，オブジェクトの干渉チェックの実装を容易にしている．

`search`前に，`this->networks`の全`networkLine`のstatusを反転させておく．
また，**this->networksにないnetworkには`searcher`は入れない．**
search_detail*/

  void search(const bool getstartObj = true) {
#ifdef debug_Network
    Print("start search", lRed);
#endif
    initSearch();
    this->depth = 0;
    this->netObjsAtDepth.clear();
    this->netObjsAtDepthUO.clear();

    // for (const auto &obj : this->getObjects())
    // 	obj->remember(this);

    V_netPp nextP = netObjs /*{this->startObj}*/, checkP;
    this->netObjsAtDepth.emplace_back(nextP);

    this->netObjsAtDepthUO.emplace_back(netObjsUO);

    if (!getstartObj) {
      this->netObjs_.clear();
      auto temp = this->netObjs;
      this->netObjs = this->netObjs_;
      this->netObjs_ = temp;
      /* ------------------------------------------------------ */
      this->netObjs_UO.clear(); // swapを使おう
      auto tempUO = this->netObjsUO;
      this->netObjsUO = this->netObjs_UO;
      this->netObjs_UO = tempUO;
    }

    networkPoint *P;
    do {
      checkP = nextP;
      nextP.clear();

      this->depth++;
      // V_netPp netObjsThisDepth = {};
      this->netObjsAtDepth.push_back({});
      this->netObjsAtDepthUO.push_back({});
      for (const auto &p : checkP) {
        for (const auto &l : p->getLines()) {
          if (!this->enteredLines_UO.count(l) && !this->enteredLinesUO.count(l) /*!l->doUKM(this)*/ /*mandetory*/) {
            if (this->condEnterLine(p, l)) {

              this->enteredLines.push_back(l);
              this->enteredLinesUO.emplace(l);
              if ((P = (*l)(p)) != nullptr /*mandetory*/ && !this->netObjsUO.count(P) && !this->netObjs_UO.count(P) /*!P->doUKM(this)*/ /*mandetory*/) {
                this->reachedLines.push_back(l);
                this->reachedLinesUO.emplace(l);
                if (this->condGetObject(l, P)) {
                  // netObjsThisDepth.emplace_back(P);
                  this->netObjsAtDepth.rbegin()->emplace_back(P); // 最後尾のベクトルに追加していく
                  this->netObjs.push_back(P);
                  /* ------------------------------------------------------ */
                  this->netObjsAtDepthUO.rbegin()->emplace(P); // 最後尾のベクトルに追加していく
                  this->netObjsUO.emplace(P);
                } else {
                  this->netObjs_.push_back(P);
                  this->netObjs_UO.emplace(P);
                }
                if (this->condKeepSearch(l, P)) {
                  nextP.emplace_back(P);
                }
              } else {
                this->reachedLines_.push_back(l);
                this->reachedLines_UO.emplace(l);
              }
            } else {
              this->enteredLines_.push_back(l);
              this->enteredLines_UO.emplace(l);
            }
          }
        }
      }
      // this->netObjsAtDepth.emplace_back(netObjsThisDepth);
    } while (!nextP.empty());

    std::sort(this->enteredLines.begin(), this->enteredLines.end());
    std::sort(this->enteredLines_.begin(), this->enteredLines_.end());
    std::sort(this->reachedLines.begin(), this->reachedLines.end());
    std::sort(this->reachedLines_.begin(), this->reachedLines_.end());
    std::sort(this->netObjs.begin(), this->netObjs.end());
    std::sort(this->netObjs_.begin(), this->netObjs_.end());
  };
};
////////////////////////////////////////////////
template <> class searcher<networkFace> : public searcherCommon<networkFace> {
public:
  std::vector<T2Tddd> path;
  searcher() : searcherCommon<networkFace>() {};

  /*search_detail
点と点は，直接繋がっておらず，間には`networkLine`が存在する．これによって，`search()`は，`condEnterLine`を使って，通過する線の条件も変更することができる．このことは，オブジェクトの干渉チェックの実装を容易にしている．

`search`前に，`this->networks`の全`networkLine`のstatusを反転させておく．
また，**this->networksにないnetworkには`searcher`は入れない．**
search_detail*/
  // #define debug_searcher_networkFace
  void search(const bool getstartObj = true) {
    initSearch();
    this->depth = 0;
    this->netObjsAtDepth.clear();
    this->netObjsAtDepthUO.clear();

    // for (const auto &obj : this->getObjects())
    // 	obj->remember(this);

    std::vector<networkFace *> nextP = netObjs /*{this->startObj}*/, checkP;
    this->netObjsAtDepth.emplace_back(nextP);

    this->netObjsAtDepthUO.emplace_back(netObjsUO);

    if (!getstartObj) {
      this->netObjs_.clear();
      auto temp = this->netObjs;
      this->netObjs = this->netObjs_;
      this->netObjs_ = temp;
      /* ------------------------------------------------------ */
      this->netObjs_UO.clear(); // swapを使おう
      auto tempUO = this->netObjsUO;
      this->netObjsUO = this->netObjs_UO;
      this->netObjs_UO = tempUO;
    }

    // networkFace *P;
    V_netFp Ps;
    do {
      checkP = nextP;
      nextP.clear();
      //------------
      this->depth++;
      std::vector<networkFace *> netObjsThisDepth = {};
      std::unordered_set<networkFace *> netObjsThisDepthUO = {};
      for (const auto &p : checkP) {
        // for (const auto &l : p->Lines)
        std::ranges::for_each(p->Lines, [&](const auto &l) {
          if (!this->enteredLines_UO.count(l) && !this->enteredLinesUO.count(l) /*!l->doUKM(this)*/ /*mandetory*/) {
            if (this->condEnterLine(p, l)) {
              // this->enteredLines.emplace_back(l);

              this->enteredLines.push_back(l);
              this->enteredLinesUO.emplace(l);
              Ps = l->getFaces();
              network::erase(Ps, p);
              for (const auto &P : Ps) {
                // 一対一の関係ではない
                if ((P != nullptr) /*mandetory*/ && !this->netObjsUO.count(P) && !this->netObjs_UO.count(P) /*!P->doUKM(this)*/ /*mandetory*/) {
                  this->reachedLines.push_back(l);
                  this->reachedLinesUO.emplace(l);
                  if (this->condGetObject(l, P)) {
                    // this->netObjs.emplace_back(P);
                    this->netObjs.push_back(P);
                    netObjsThisDepth.emplace_back(P);
                    path.push_back({Mean(p->getLocationsTuple()), Mean(P->getLocationsTuple())});
                    /* ------------------------------------------------------ */
                    this->netObjsUO.emplace(P);
                    netObjsThisDepthUO.emplace(P);
                  } else {
                    this->netObjs_.push_back(P);
                    this->netObjs_UO.emplace(P);
                  }
                  if (this->condKeepSearch(l, P)) {
                    nextP.emplace_back(P);
                  }
                } else {
                  // this->reachedLines_.emplace_back(l);
                  this->reachedLines_.push_back(l);
                  this->reachedLines_UO.emplace(l);
                }
              }
            } else {
              // this->enteredLines_.emplace_back(l);
              this->enteredLines_.push_back(l);
              this->enteredLines_UO.emplace(l);
            }
          };
        });
      }
      this->netObjsAtDepth.emplace_back(netObjsThisDepth);
      this->netObjsAtDepthUO.emplace_back(netObjsThisDepthUO);
    } while (!nextP.empty());

    std::sort(this->enteredLines.begin(), this->enteredLines.end());
    std::sort(this->enteredLines_.begin(), this->enteredLines_.end());
    std::sort(this->reachedLines.begin(), this->reachedLines.end());
    std::sort(this->reachedLines_.begin(), this->reachedLines_.end());
    std::sort(this->netObjs.begin(), this->netObjs.end());
    std::sort(this->netObjs_.begin(), this->netObjs_.end());
  };
};

//-----------------------------------------------------------
template <typename T> class depth_searcher : public searcher<T> {
public:
  int depth_lim;
  depth_searcher(const int depth_lim_IN) : depth_lim(depth_lim_IN), searcher<T>() {};
  depth_searcher() : searcher<T>() {};
  void setDepth(const int depth_lim_IN) { this->depth_lim = depth_lim_IN; };
  bool condEnterLine(const T *p, const netLp l) override {
    // 自身(startobj):0
    //  neigbors:1
    if (this->depth <= depth_lim /*0なら検索なし*/)
      return true;
    else
      return false;
  };
  bool condGetObject(const netLp l, const T *P) override {
    if (MemberQ(this->networks, P->getNetwork()))
      return true;
    else
      return false;
  };
};
//-----------------------------------------------------------
// 名前がしっくりこなかったので新たな名前で作り直した
template <typename T> class BreadthFirstSearcher : public searcher<T> {
protected:
  int depth_lim;

public:
  BreadthFirstSearcher(const int depth_lim_IN /*0は自分だけ，depth_limまで検索*/) : depth_lim(depth_lim_IN), searcher<T>() {};

  std::vector<T *> operator()(const int n) {
    if (n >= static_cast<int>(this->netObjsAtDepth.size())) {
      std::stringstream ss;
      ss << "n = " << n << std::endl;
      ss << "netObjsAtDepth.size() = " << this->netObjsAtDepth.size() << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
    }
    return this->netObjsAtDepth[n];
  };

  std::vector<T *> operator()() { return this->netObjsAtDepth; };

  // void setDepth(const int depth_lim_IN) { this->depth_lim = depth_lim_IN; };
  bool condEnterLine(const T *p, const netLp l) override {
    if (this->depth <= depth_lim /*0なら検索なし．defaultで自分を含める．初めのループではdepth=1*/)
      return true;
    else
      return false;
  };
  bool condGetObject(const netLp l, const T *P) override {
    if (MemberQ(this->networks, P->getNetwork()))
      return true;
    else
      return false;
  };
  bool condKeepSearch(const netLp l, const T *P) override {
    if (this->depth == depth_lim /*0なら検索なし．defaultで自分を含める．初めのループではdepth=1*/)
      return false;
    else
      return true;
  };
};
//////////////////////////////////////////
inline VV_netPp BFS(networkPoint *startObj, const int depth, const V_Netp &nets = {}) {
  BreadthFirstSearcher<networkPoint> bfs(depth);
  bfs.set(startObj);
  bfs.addNetworks(nets);
  bfs.search();
  return bfs.getObjectsAtDepth();
};
inline VV_netPp BFS(const V_netPp &startObj, const int depth, const V_Netp &nets = {}) {
  BreadthFirstSearcher<networkPoint> bfs(depth);
  bfs.set(startObj);
  bfs.addNetworks(nets);
  bfs.search();
  return bfs.getObjectsAtDepth();
};
inline VV_netFp BFS(const netFp startObj, const int depth, const V_Netp &nets = {}) {
  BreadthFirstSearcher<networkFace> bfs(depth);
  bfs.set(startObj);
  bfs.addNetworks(nets);
  bfs.search();
  return bfs.getObjectsAtDepth();
};
inline VV_netFp BFS(const V_netFp &startObj, const int depth, const V_Netp &nets = {}) {
  BreadthFirstSearcher<networkFace> bfs(depth);
  bfs.set(startObj);
  bfs.addNetworks(nets);
  bfs.search();
  return bfs.getObjectsAtDepth();
};
inline VV_netFp BFS(const std::unordered_set<networkFace *> &startObj, const int depth, const V_Netp &nets = {}) {
  BreadthFirstSearcher<networkFace> bfs(depth);
  bfs.set(ToVector(startObj));
  bfs.addNetworks(nets);
  bfs.search();
  return bfs.getObjectsAtDepth();
};
inline std::unordered_set<networkFace *> BFS_Flattened(const std::unordered_set<networkFace *> &startObj, const int depth, const V_Netp &nets = {}) {
  BreadthFirstSearcher<networkFace> bfs(depth);
  bfs.set(ToVector(startObj));
  bfs.addNetworks(nets);
  bfs.search();
  return bfs.getObjectsUO();
};
/* ------------------------------------------------------ */
inline std::vector<std::unordered_set<networkPoint *>> BFSUO(const netPp startObj, const int depth, const V_Netp &nets = {}) {
  BreadthFirstSearcher<networkPoint> bfs(depth);
  bfs.set(startObj);
  bfs.addNetworks(nets);
  bfs.search();
  return bfs.getObjectsAtDepthUO();
};
inline std::unordered_set<networkPoint *> BFS_Flattened(const netPp startObj, const int depth, const V_Netp &nets = {}) {
  BreadthFirstSearcher<networkPoint> bfs(depth);
  bfs.set(startObj);
  bfs.addNetworks(nets);
  bfs.search();
  return bfs.getObjectsUO();
};
inline std::vector<std::unordered_set<networkPoint *>> BFSUO(const V_netPp &startObj, const int depth, const V_Netp &nets = {}) {
  BreadthFirstSearcher<networkPoint> bfs(depth);
  bfs.set(startObj);
  bfs.addNetworks(nets);
  bfs.search();
  return bfs.getObjectsAtDepthUO();
};
inline std::vector<std::unordered_set<networkFace *>> BFSUO(const netFp startObj, const int depth, const V_Netp &nets = {}) {
  BreadthFirstSearcher<networkFace> bfs(depth);
  bfs.set(startObj);
  bfs.addNetworks(nets);
  bfs.search();
  return bfs.getObjectsAtDepthUO();
};
inline std::vector<std::unordered_set<networkFace *>> BFSUO(const V_netFp &startObj, const int depth, const V_Netp &nets = {}) {
  BreadthFirstSearcher<networkFace> bfs(depth);
  bfs.set(startObj);
  bfs.addNetworks(nets);
  bfs.search();
  return bfs.getObjectsAtDepthUO();
};

/* -------------------------------------------------------------------------- */
inline std::unordered_set<networkFace *> bfs(const std::unordered_set<networkFace *> &FACES, const short unsigned int s) {
  // if s=0 ->  do nothing
  std::unordered_set<networkFace *> tmp = FACES, ret = FACES;
  for (auto i = 0; i < s; i++) {
    for (const auto &F : tmp) {
      std::ranges::for_each(F->getPoints(), [&](const auto &p) {
        for (const auto &f : p->getFaces())
          ret.emplace(f);
      });
    }
    tmp = ret;
  }
  return ret;

  // 段階をわければ，全てを改めて使ってfor_eachしなくていいので，早く終わるのでは？

  // std::unordered_set<networkFace *> accum = FACES;
  // std::vector<networkFace *> newfaces, loopwith(FACES.begin(), FACES.end());
  // for (auto i = 0; i < s; i++) {
  //    newfaces.clear();
  //    for (const auto &F : loopwith) {
  //       for_each(F->getPoints(), [&](const auto &p) {
  //          for (const auto &f : p->getFaces()) {
  //             if (*((accum.emplace(f)).first))
  //                newfaces.emplace_back(f);  // 新たに追加された
  //          }
  //       });
  //    }
  //    loopwith = newfaces;
  // }
  // return accum;
};

inline std::unordered_set<networkFace *> bfs(const std::vector<networkFace *> &FACES, const short unsigned int s) {
  // if s=0 ->  do nothing
  std::unordered_set<networkFace *> tmp(FACES.begin(), FACES.end()), ret(FACES.begin(), FACES.end());
  for (auto i = 0; i < s; i++) {
    for (const auto &F : tmp) {
      std::ranges::for_each(F->getPoints(), [&](const auto &p) {
        for (const auto &f : p->getFaces())
          ret.emplace(f);
      });
    }
    tmp = ret;
  }
  return ret;

  // 段階をわければ，全てを改めて使ってfor_eachしなくていいので，早く終わるのでは？

  // std::unordered_set<networkFace *> accum = FACES;
  // std::vector<networkFace *> newfaces, loopwith(FACES.begin(), FACES.end());
  // for (auto i = 0; i < s; i++) {
  //    newfaces.clear();
  //    for (const auto &F : loopwith) {
  //       for_each(F->getPoints(), [&](const auto &p) {
  //          for (const auto &f : p->getFaces()) {
  //             if (*((accum.emplace(f)).first))
  //                newfaces.emplace_back(f);  // 新たに追加された
  //          }
  //       });
  //    }
  //    loopwith = newfaces;
  // }
  // return accum;
};

#endif
