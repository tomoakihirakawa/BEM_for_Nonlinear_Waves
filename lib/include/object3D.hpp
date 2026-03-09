#ifndef INCL_object3D
#define INCL_object3D

using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;
using V_int = std::vector<int>;
using VV_int = std::vector<std::vector<int>>;

//=====================================================================
//============================= object3D ==============================
//=====================================================================
namespace obj3D {
using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;
//-----------------------------------

// 境界面の重なりを検知するのみ，完全に内部に入った場合はfalse
bool isBoundariesOverlap(const VV_d &boundsA, const VV_d &boundsB) {
   // for (auto i = 0; i < 3 /*rboundsまで広げるとsetCrossPointsがうまくいかない*/; i++)
   // {
   //   // minA   maxA    minB    maxB
   //   //  *-------*  <  *-------*
   //   // minB   maxB  minA    maxA
   //   //  *-------*  <  *-------*
   //   if (boundsA[i][1] /*max*/ < boundsB[i][0] /*min*/ || boundsB[i][1] /*max*/ < boundsA[i][0] /*min*/)
   //     return false;
   // }
   if ((boundsA[0][1] /*max*/ < boundsB[0][0] /*min*/ || boundsB[0][1] /*max*/ < boundsA[0][0] /*min*/) ||
       (boundsA[1][1] /*max*/ < boundsB[1][0] /*min*/ || boundsB[1][1] /*max*/ < boundsA[1][0] /*min*/) ||
       (boundsA[2][1] /*max*/ < boundsB[2][0] /*min*/ || boundsB[2][1] /*max*/ < boundsA[2][0] /*min*/))
      return false;

   return true;
};
//--- extrat ---
// vector内のオブジェクトの内部の情報をvectorとして返す
bool isBoundaryInside(const VV_d &boundsA, const VV_d &boundsB) {
   for (auto i = 0; i < 3 /*rboundsまで広げるとsetXPointsがうまくいかない*/; i++) {
      // minA   minB    maxB  maxA
      //  *      *-------*     *
      if (!(boundsA[i][0] /*min*/ <= boundsB[i][0] /*min*/ && boundsB[i][1] /*min*/ <= boundsA[i][1] /*max*/))
         return false;  // B is not inside
   }
   return true;
};
//--------------------------------
template <class T>
VV_d extractX(const std::vector<T *> &object) {
   VV_d ret(object.size(), V_d(3));
   for (auto i = 0; i < object.size(); i++)
      ret[i] = object[i]->getX();
   return ret;
};
template <class T>
VVV_d extractX(const std::vector<std::vector<T *>> &object) {
   VVV_d ret;
   ret.reserve(object.size());
   for (const auto &obj : object)
      ret.emplace_back(obj3D::extractX(obj));
   return ret;
};
//--------------------------------
template <class T>
std::unordered_set<T *> takeIfBoundariesOverlap(const std::unordered_set<T *> &objs, const CoordinateBounds &bounds) {
   std::unordered_set<T *> ret;
   for (const auto &obj : objs)
      if (IntersectQ(obj->getBounds(), bounds))
         ret.emplace(obj);
   return ret;
};
template <class T>
std::vector<T *> takeIfBoundariesOverlap(const std::vector<T *> &objs, const CoordinateBounds &bounds) {
   std::vector<T *> ret;
   for (const auto &obj : objs)
      if (IntersectQ(obj->getBounds(), bounds))
         ret.emplace_back(obj);
   return ret;
};
template <class T>
std::vector<T *> takeIfBoundariesOverlap(const std::vector<T *> &objs, const VV_d &bounds) {
   std::vector<T *> ret;
   for (const auto &obj : objs)
      if (obj3D::isBoundariesOverlap(obj->getBounds(), bounds))
         ret.emplace_back(obj);
   return ret;
};
template <class T>
std::vector<T *> takeIfNotBoundariesOverlap(const std::vector<T *> &objs, const VV_d &bounds) {
   std::vector<T *> ret;
   for (const auto &obj : objs)
      if (!obj3D::isBoundariesOverlap(obj->getBounds(), bounds))
         ret.emplace_back(obj);
   return ret;
};
template <class T>
std::vector<T *> takeInsideOfBounds(const std::vector<T *> &objs, const VV_d &bounds) {
   std::vector<T *> ret;
   for (const auto &obj : objs)
      if (obj3D::isBoundaryInside(bounds, obj->getBounds()))
         ret.emplace_back(obj);
   return ret;
};
}  // namespace obj3D

struct object3D {
   object3D(const CoordinateBounds &bs) {
      this->bounds = bs;
      std::get<0>(this->X) = Mean(std::get<0>(bs.bounds));
      std::get<1>(this->X) = Mean(std::get<1>(bs.bounds));
      std::get<2>(this->X) = Mean(std::get<2>(bs.bounds));
      if (!isFinite(X)) {
         std::stringstream ss;
         ss << "this->X = " << this->X << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      };
   };
   ~object3D(){};
   /* -------------------------------------------- */
   Tddd X;  // 新たにxyzの代わりに使っていく座標2021/08/23
   const Tddd &getXtuple() const { return this->X; /*面のsetBoundsでは，T3Tdddの平均がXとなるようになっている．バウンディングボックスの中心ではない.*/ };
   CoordinateBounds bounds;
   //
   V_d getX() const { return {std::get<0>(this->X), std::get<1>(this->X), std::get<2>(this->X)}; };
   const CoordinateBounds &getBounds() const { return this->bounds; };
   double getScale() const {
      return Norm(Tddd{std::get<1>(std::get<0>(this->bounds.bounds)) - std::get<0>(std::get<0>(this->bounds.bounds)),
                       std::get<1>(std::get<1>(this->bounds.bounds)) - std::get<0>(std::get<1>(this->bounds.bounds)),
                       std::get<1>(std::get<2>(this->bounds.bounds)) - std::get<0>(std::get<2>(this->bounds.bounds))});
   };
   void setBounds(const CoordinateBounds &bs) {
      this->bounds = bs;
      std::get<0>(this->X) = Mean(std::get<0>(bs.bounds));
      std::get<1>(this->X) = Mean(std::get<1>(bs.bounds));
      std::get<2>(this->X) = Mean(std::get<2>(bs.bounds));
   };
   void setBounds(const T3Tddd &X012) {
      this->setBounds(CoordinateBounds(X012));
      this->X = Mean(X012);
   };
};
#endif
