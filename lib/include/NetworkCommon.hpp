#ifndef NetworkCommon_H
#define NetworkCommon_H

#include "basic.hpp"
#include <cstdio>
#include <map>
#include <vector>
namespace network { // tempalte fucntions which are useed for network object vectors

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;
using V_int = std::vector<int>;
using VV_int = std::vector<std::vector<int>>;

template <typename T> int find(const std::vector<T *> &P, const T *p) {
  for (int i = 0; i < P.size(); i++)
    if (P[i] == p)
      return i;
  return -1;
};
//
template <typename T> bool add(std::vector<T *> &P, T *p) {
  if (P.empty()) {
    P.emplace_back(p);
    return true;
  }
  auto it = std::find(P.begin(), P.end(), p);
  if (it == P.end()) {
    P.emplace_back(p);
    return true;
  }
  return false;
};
//
template <typename T> std::vector<bool> add(std::vector<T *> &P, const std::vector<T *> &p) {
  std::vector<bool> ret;
  ret.reserve(p.size());
  for (const auto &q : p)
    ret.emplace_back(network::add(P, q));
  return ret;
};
//-------------------------------
template <typename T> bool erase(std::vector<T *> &P, const T *p) {
  auto it = std::find(P.begin(), P.end(), p);
  if (it != P.end()) {
    P.erase(it);
    return true;
  }
  return false;
};
//-------------------------------
// TakeExceptで置き換える
// template <typename T>
// std::vector<bool> erase(std::vector<T *> &P, const std::vector<T *> &p)
// {
//   std::vector<bool> ret;
//   ret.reserve(p.size());
//   for (const auto &q : p)
//     ret.emplace_back(network::erase(P, q));
//   return ret;
// };
//-------------------------------
template <typename T, typename U> bool erase(std::map<T *, U> &P, T *p) {
  auto it = P.find(p);
  if (it != P.end()) {
    P.erase(it);
    return true;
  }
  return false;
};
//////////////////////////////////
//-------------------------------
template <typename T> void void_erase(std::vector<T *> &P, const T *p) {
  auto it = std::find(P.cbegin(), P.cend(), p);
  if (it != P.end())
    P.erase(it);
};
//-------------------------------
template <typename T> void void_erase(std::vector<T *> &P, const std::vector<T *> &ps) {
  for (const auto &q : ps)
    network::void_erase(P, q);
};
//-------------------------------
template <typename T> bool myswitch(std::vector<T *> &L, const T *oldL, T *newL) {
  //   注意：重複を許す事になっている
  //  +----+    +----+    +----+    +----+
  //  |newL|    |  F-|--->|oldL|    |newF|
  //  |    |    |    |<---|-   |    |    |
  //  +----+    +----+    +----+    +----+
  for (auto &l : L)
    if (l == oldL) {
      l = newL;
      //  +----+    +----+    +----+    +----+
      //  |newL|<---|- F |    |oldL|    |newF|
      //  |    |    |    |<---|-   |    |    |
      //  +----+    +----+    +----+    +----+
      return true;
    }
  return false;
};

/*
replaceは，Fに接続しているoldLをFから剥ぎ取り，newObjに付け替える．
oldObjはFの代わりにnewFに接続しなおす．
*/
// template <typename L, typename T>
// bool replace(T *this_, std::vector<L *> &this_V_L, L *oldL, L *newL, T *newF = nullptr)
// { /*this_ may be holded by L*/
//   //  +----+    +----+    +----+    +----+
//   //  |newL|    |this|--->|oldL|    |newF|
//   //  |    |    |PorF|<---|-   |    |    |
//   //  +----+    +----+    +----+    +----+
//   if (myswitch(this_V_L, oldL, newL) /*if found*/)
//   {
//     //  +----+    +----+    +----+    +----+
//     //  |newL|<---|this|    |oldL|    |newF|
//     //  |    |    |PorF|<---|-   |    |    |
//     //  +----+    +----+    +----+    +----+
//     oldL->erase(this_ /*PorF*/);
//     //  +----+    +----+    +----+    +----+
//     //  |newL|<---|this|    |oldL|    |newF|
//     //  |    |    |PorF|    |    |    |    |
//     //  +----+    +----+    +----+    +----+
//     newL->Add(this_ /*PorF*/);
//     //  +----+    +----+    +----+    +----+
//     //  |newL|<---|this|    |oldL|    |newF|
//     //  |   -|--->|PorF|    |    |    |    |
//     //  +----+    +----+    +----+    +----+
//     if (newF != nullptr)
//     {
//       oldL->Add(newF /*PorF*/);
//       //  +----+    +----+    +----+    +----+
//       //  |newL|<---|this|    |oldL|    |newF|
//       //  |   -|--->|PorF|    |   -|--->|    |
//       //  +----+    +----+    +----+    +----+
//       newF->Add(oldL /*PorF*/);
//       //  +----+    +----+    +----+    +----+
//       //  |newL|<---|this|    |oldL|<---|newF|
//       //  |   -|--->|PorF|    |   -|--->|    |
//       //  +----+    +----+    +----+    +----+
//     }
//     return true;
//   }
//   return false;
// };

//////////////////////////////////////////
// template <typename T, typename U>
// bool doubleReplace(U *Uthis, T *Ta /*孤立するオブジェクト*/, T *Tb, U *Uother)
// {
//   Uthis->swap(Ta, Tb); //1
//   Ta->erase(Uthis);      //2
//   Tb->Add(Uthis);        //3
//                          //このステップがdouble replace
//                          //switchでないと，順番に意味のあるFaceではおかしくなるので注意
//   Ta->Add(Uother);
//   Uother->Add(Ta);
//   // throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//   return true;
// };
/////////////////////////////////////
template <class T> std::vector<T *> sortByDistance(const std::vector<T *> &Ps, T *P) {
  std::vector<T *> ret;
  bool inseted;
  for (const auto &p : Ps) {
    inseted = false;
    for (auto it = ret.begin(); it != ret.end(); ++it) {
      if (Norm(P->getXtuple() - p->getXtuple()) < Norm(P->getXtuple() - (*it)->getXtuple())) {
        ret.insert(it, p);
        inseted = true;
        break;
      }
    }
    if (!inseted)
      ret.emplace_back(p);
  }
  return ret;
};

template <class T> std::vector<T *> ordering(std::vector<T *> remain) {
  //    V_netPp remain = getXPointsIntersect();
  std::vector<T *> ret = {*remain.begin()};
  remain.erase(remain.begin());
  bool found;
  while (remain.size() > 0) {
    found = false;
    for (auto it = remain.begin(); it != remain.end(); it++) {
      if (network::find((*ret.begin())->getNeighbors(), *it) != -1) {
        ret.insert(ret.begin(), *it);
        remain.erase(it);
        found = true;
        break;
      } else if (network::find((*ret.rbegin())->getNeighbors(), *it) != -1) {
        ret.insert(ret.end(), *it);
        remain.erase(it);
        found = true;
        break;
      }
    }
    if (!found) {
      Print(__PRETTY_FUNCTION__, Red);
      Print("not found", red);
      Print(remain, Red);
      Print(ret, blue);
      ret.insert(ret.end(), remain.begin(), remain.end());
      std::cin.ignore();
      break;
    }
  }
  return ret;
};
//----------------------
//------- status -------
//----------------------
template <class T> void setStatus(const std::vector<T *> &Objects, bool TorF) {
  for (const auto &object : Objects)
    object->setStatus(TorF);
};
template <class T> void setStatus(const std::unordered_set<T *> &Objects, bool TorF) {
  for (const auto &object : Objects)
    object->setStatus(TorF);
};
template <class T> std::vector<T *> getIfStatus(const std::vector<T *> &Objects, bool TorF = true) {
  std::vector<T *> ret;
  for (const auto &object : Objects)
    if (object->getStatus() == TorF)
      ret.emplace_back(object);
  return ret;
};

template <class T> bool AnyStatusTrue(const std::vector<T *> &obj) {
  for (const auto &o : obj)
    if (o->getStatus())
      return true;
  return false;
};

template <class T> bool AllStatusTrue(const std::vector<T *> &obj) {
  for (const auto &o : obj)
    if (!o->getStatus())
      return false;
  return true;
};
//-----------------------
//-------- take ---------
//-----------------------
// vectorの中から特定のものを抜き出す関数
template <class T> T *takeNearest(const std::vector<T *> &objs, const Tddd &xyz) {
  if (objs.empty()) {
    Print(__PRETTY_FUNCTION__, Red);
    abort();
  }
  T *ret;
  double min = 1E+30, dist;
  for (const auto &obj : objs) {
    dist = Norm(obj->getXtuple() - xyz);
    if (min > dist) {
      ret = obj;
      min = dist;
    }
  }
  return ret;
};
template <class T> T *takeNearest(const std::unordered_set<T *> &objs, const Tddd &xyz) {
  if (objs.empty()) {
    Print(__PRETTY_FUNCTION__, Red);
    abort();
  }
  T *ret;
  double min = 1E+30, dist;
  for (const auto &obj : objs) {
    dist = Norm(obj->getXtuple() - xyz);
    if (min > dist) {
      ret = obj;
      min = dist;
    }
  }
  return ret;
};
template <class T> std::vector<T *> takeIfStatus(const std::vector<T *> &obj, const bool TorF = true) {
  std::vector<T *> ret;
  for (const auto &o : obj) {
    if (o->getStatus() == TorF) {
      ret.emplace_back(o);
    }
  }
  return ret;
};

template <class T> std::vector<T *> takeIfIntersect(const std::vector<T *> &obj, const bool TorF = true) {
  std::vector<T *> ret;
  for (const auto &f : obj) {
    if (f->intersectQ() == TorF) {
      ret.emplace_back(f);
    }
  }
  return ret;
};

template <class T, class Network> std::vector<T *> takeIfNetwork(const std::vector<T *> &obj, const Network *net) {
  std::vector<T *> ret;
  for (const auto &f : obj) {
    if (f->getNetwork() == net) {
      ret.emplace_back(f);
    }
  }
  return ret;
};

} // namespace network

///////////////////////////////////////////////////////////
template <typename T> bool isEdge(const T *l) { // 周囲の全ての線が面を2つ保持している場合，false
  return (l->getFaces().size() < 2);
};

inline bool isEdgePoint(const auto &p) {
  return p->getLines().empty() || std::any_of(p->getLines().begin(), p->getLines().end(), [&](const auto &l) { return (!l || l->getFaces().size() <= 1); });
};

#endif
