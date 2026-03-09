#ifndef lib_computationalGeometry_H
#define lib_computationalGeometry_H
#include "basic_vectors.hpp"
/*
circular dependenciesを認めよう．
ただし，circular dependenctなクラス自体が大きくなると，
保守が大変になるので（一つの修正が全てに及ぶため），
できるだけ，派生させた上，末端に作成する方が良いと思う．
また，メンバー関数を減らし，外部に関数を作る方がいい．
*/

namespace geometric_network
{
    /* ------------------------------------------------------ */
    class Point;
    class Edge;
    class Triangle;
    class Tetrahedron;
    /* ------------------------------------------------------ */
    class Network
    {
    public:
        std::unordered_set<Point *> points;
        std::unordered_set<Edge *> edges;
        std::unordered_set<Triangle *> triangles;
        std::unordered_set<Tetrahedron *> tetras;
        void add(Point *p) { this->points.emplace(p); };
        void add(Edge *p) { this->edges.emplace(p); };
        void add(Triangle *p) { this->triangles.emplace(p); };
        void add(Tetrahedron *p) { this->tetras.emplace(p); };
        Network(){};
        Edge *findEdge(Point *p0, Point *p1);
        Triangle *findTriangle(Point *p0, Point *p1, Point *p2);
    };
    /* ------------------------------------------------------ */
    class Point
    {
    public:
        Tddd X;
        Network *network;
        std::unordered_set<Edge *> edges;
        std::unordered_set<Triangle *> triangles;
        std::unordered_set<Tetrahedron *> tetras;
        void add(Edge *p) { this->edges.emplace(p); };
        void add(Triangle *p) { this->triangles.emplace(p); };
        void add(Tetrahedron *p) { this->tetras.emplace(p); };
        //
        Point(Network *net, const Tddd &XIN) : network(net), X(XIN) { net->add(this); };
        Point(Network *net, const double &x, const double &y, const double &z) : X({x, y, z}) { net->add(this); };
    };
    /* ------------------------------------------------------ */
    class Edge
    {
    public:
        Network *network;
        Tddd X;
        Point *P0, *P1;
        std::unordered_set<Triangle *> triangles;
        std::unordered_set<Tetrahedron *> tetras;
        void add(Triangle *p) { this->triangles.emplace(p); };
        void add(Tetrahedron *p) { this->tetras.emplace(p); };
        Edge(Network *net, Point *p0, Point *p1) : network(net), P0(p0), P1(p1), X(0.5 * (p0->X + p1->X))
        {
            std::cout << p0 << ", " << p1 << std::endl;
            P0->edges.emplace(this);
            P1->edges.emplace(this);
            net->add(this);
        };
        std::tuple<Point *, Point *> getPoints() const { return {P0, P1}; };
    };
    /* ------------------------------------------------------ */
    std::set<Point *> sortedPoints(Edge *E0, Edge *E1, Edge *E2)
    {
        std::set<Point *> ret;
        if (E2->P0 == E0->P0 || E2->P1 == E0->P0)
        {
            ret.insert(E0->P0);
            ret.insert(E0->P1);
        }
        else if (E2->P0 == E0->P1 || E2->P1 == E0->P1)
        {
            ret.insert(E0->P1);
            ret.insert(E0->P0);
        }
        ret.insert(E1->P0);
        ret.insert(E1->P1);
        ret.insert(E2->P0);
        ret.insert(E2->P1);
        return ret;
        // usage:
        // {*it, *std::next(it, 1), *std::next(it, 2)};
    };
    /* ------------------------------------------------------ */
    class Triangle
    {
    public:
        Network *network;
        Tddd X;
        Point *P0, *P1, *P2; // ordered
        Edge *E0, *E1, *E2;  // ordered
        Tetrahedron *Tetra0, *Tetra1;
        Tetrahedron *add(Tetrahedron *t)
        {
            //この三角形から取り除かれるテトラを返す
            auto ret = Tetra1;
            Tetra1 = Tetra0;
            Tetra0 = t;
            return ret;
        };
        //
        std::tuple<Point *, Point *, Point *> getPoints() const { return {P0, P1, P2}; };
        std::tuple<Edge *, Edge *, Edge *> getEdges() const { return {E0, E1, E2}; };
        std::tuple<Tetrahedron *, Tetrahedron *> getTetras() const { return {Tetra0, Tetra1}; };
        //
        Triangle(Network *net, Point *p0, Point *p1, Point *p2)
            : network(net),
              X((p0->X + p1->X + p2->X) / 3.),
              P0(p0), P1(p1), P2(p2),
              E0(nullptr), E1(nullptr), E2(nullptr),
              Tetra0(nullptr), Tetra1(nullptr)
        {
            std::cout << P0 << ", " << P1 << ", " << P2 << std::endl;
            P0->add(this);
            P1->add(this);
            P2->add(this);
            net->add(this);
        };
        Triangle(Network *net, Edge *e0, Edge *e1, Edge *e2)
            : network(net),
              X((e0->X + e1->X + e2->X) / 3.),
              P0(nullptr), P1(nullptr), P2(nullptr),
              E0(e0), E1(e1), E2(e2),
              Tetra0(nullptr), Tetra1(nullptr)
        {
            E0->add(this);
            E1->add(this);
            E2->add(this);
            auto sp = sortedPoints(E0, E1, E2);
            auto it = sp.begin();
            this->P0 = *it;
            this->P1 = *std::next(it, 1);
            this->P2 = *std::next(it, 2);
            net->add(this);
        };
    };
    /* ------------------------------------------------------ */
    class Tetrahedron
    {
        //     1,3,2
        //          3
        // 0,2,3  / | \ 0,3,1     --- 1 ---
        //       2--|--1           \      /
        // 0,1,2  \ | /             2    0
        //          0                 \/
    public:
        Network *network;
        Tddd X;
        Point *P0, *P1, *P2, *P3;
        Edge *E0, *E1, *E2, *E3, *E4, *E5; // ordered
        Triangle *T0, *T1, *T2, *T3;
        //
        Tetrahedron(Network *net)
            : network(net),
              P0(nullptr), P1(nullptr), P2(nullptr), P3(nullptr),
              E0(nullptr), E1(nullptr), E2(nullptr), E3(nullptr), E4(nullptr), E5(nullptr),
              T0(nullptr), T1(nullptr), T2(nullptr), T3(nullptr)
        {
            std::cout << P0 << ", " << P1 << ", " << P2 << ", " << P3 << std::endl;
            net->add(this);
        };
        Tetrahedron(Network *net, Point *p0, Point *p1, Point *p2, Point *p3)
            : network(net),
              X((p0->X + p1->X + p2->X + p3->X) / 4.),
              P0(p0), P1(p1), P2(p2), P3(p3),
              E0(nullptr), E1(nullptr), E2(nullptr), E3(nullptr), E4(nullptr), E5(nullptr),
              T0(nullptr), T1(nullptr), T2(nullptr), T3(nullptr)
        {
            std::cout << P0 << ", " << P1 << ", " << P2 << ", " << P3 << std::endl;
            net->add(this);
            P0->add(this);
            P1->add(this);
            P2->add(this);
            P3->add(this);
        };
        //
        Tetrahedron(Network *net, Edge *e0, Edge *e1, Edge *e2, Edge *e3, Edge *e4, Edge *e5)
            : network(net),
              X((e0->X + e1->X + e2->X + e3->X + e4->X + e5->X) / 6.),
              P0(nullptr), P1(nullptr), P2(nullptr), P3(nullptr),
              E0(e0), E1(e1), E2(e2), E3(e3), E4(e4), E5(e5),
              T0(nullptr), T1(nullptr), T2(nullptr), T3(nullptr)
        {
            std::cout << E0 << ", " << E1 << ", " << E2 << ", " << E3 << ", " << E4 << ", " << E5 << std::endl;
            net->add(this);
            E0->add(this);
            E1->add(this);
            E2->add(this);
            E3->add(this);
            E4->add(this);
            E5->add(this);
        };
        //
        Tetrahedron(Network *net, Triangle *t0, Triangle *t1, Triangle *t2, Triangle *t3)
            : network(net),
              X((t0->X + t1->X + t2->X + t3->X) / 4.),
              P0(nullptr), P1(nullptr), P2(nullptr), P3(nullptr),
              E0(nullptr), E1(nullptr), E2(nullptr), E3(nullptr), E4(nullptr), E5(nullptr),
              T0(t0), T1(t1), T2(t2), T3(t3)
        {
            std::cout << P0 << ", " << P1 << ", " << P2 << ", " << P3 << std::endl;
            net->add(this);
            T0->add(this);
            T1->add(this);
            T2->add(this);
            T3->add(this);
            /*この符号で点の順番が決まる．順方向で作るcrossが外向き法線ベクトルとなるか内向きになるか*/
            if (Dot(Cross(T0->P1->X - T0->P0->X, T0->P2->X - T0->P0->X), T0->X - this->X) < 0) //テトラのP0->P1->P2で作る右手系で考える法線ベクトルは内向きになる
            {
                this->P0 = T0->P0;
                this->P1 = T0->P1;
                this->P2 = T0->P2;
            }
            else
            {
                this->P0 = T0->P2;
                this->P1 = T0->P1;
                this->P2 = T0->P0;
            }
            for_each(T1->getPoints(), [&](const auto &p)
                     {if (this->P0 != p && this->P1 != p && this->P2 != p){this->P3 = p;} });
        };
        //
        void fill(const std::tuple<Point *, Point *, Point *, Point *> &p0123)
        {
            this->fill(std::get<0>(p0123), std::get<1>(p0123), std::get<2>(p0123), std::get<3>(p0123));
        };
        void fill(Point *p0, Point *p1, Point *p2, Point *p3)
        {
            network->add(this);
            this->P0 = p0;
            this->P1 = p1;
            this->P2 = p2;
            this->P3 = p3;
            P0->add(this);
            P1->add(this);
            P2->add(this);
            P3->add(this);
            fillEdge(p0, p1, p2, p3);
            fillTriangle(p0, p1, p2, p3);
            this->X = (P0->X + P1->X + P2->X + P3->X) / 4.;
        };
        void fill(Edge *e0, Edge *e1, Edge *e2, Edge *e3, Edge *e4, Edge *e5)
        {
            auto points = sortedPoints(e0, e1, e2);
            points.insert(e3->P0);
            points.insert(e3->P1);
            auto it = points.begin();
            this->P0 = *it;
            this->P1 = *std::next(it, 1);
            this->P2 = *std::next(it, 2);
            this->P3 = *std::next(it, 3);
            fillTriangle(P0, P1, P2, P3);
        };
        void fillEdge(Point *p0, Point *p1, Point *p2, Point *p3)
        {
            // 与えられたPointを繋ぐEdgeがネットワークに存在するか確認し，あれば既存のものを利用し，なければここで生成する
            auto e = network->findEdge(p0, p1);
            E0 = (e ? e : (new Edge(network, p0, p1)));
            e = network->findEdge(p1, p2);
            E1 = (e ? e : (new Edge(network, p1, p2)));
            e = network->findEdge(p2, p0);
            E2 = (e ? e : (new Edge(network, p2, p0)));
            e = network->findEdge(p0, p3);
            E3 = (e ? e : (new Edge(network, p0, p3)));
            e = network->findEdge(p1, p3);
            E4 = (e ? e : (new Edge(network, p1, p3)));
            e = network->findEdge(p2, p3);
            E5 = (e ? e : (new Edge(network, p2, p3)));
        };
        void fillTriangle(Point *p0, Point *p1, Point *p2, Point *p3)
        {
            // 与えられたPointを繋ぐTriangleがネットワークに存在するか確認し，あれば既存のものを利用し，なければここで生成する
            auto t = network->findTriangle(p0, p1, p2);
            T0 = (t ? t : (new Triangle(network, p0, p1, p2)));
            t = network->findTriangle(p0, p3, p1);
            T1 = (t ? t : (new Triangle(network, p0, p3, p1)));
            t = network->findTriangle(p1, p2, p3);
            T2 = (t ? t : (new Triangle(network, p1, p3, p2)));
            t = network->findTriangle(p0, p2, p3);
            T3 = (t ? t : (new Triangle(network, p0, p2, p3)));
        };
        std::tuple<Point *, Point *, Point *, Point *> getPoints() const { return {P0, P1, P2, P3}; };
        std::tuple<Edge *, Edge *, Edge *, Edge *> getEdges() const { return {E0, E1, E2, E3}; };
        std::tuple<Triangle *, Triangle *, Triangle *, Triangle *> getTriangles() const { return {T0, T1, T2, T3}; };
    };
    /* ------------------------------------------------------ */
    Edge *Network::findEdge(Point *p0, Point *p1)
    {
        for (const auto &e : this->edges)
            if ((e->P0 == p0 && e->P1 == p1) ||
                (e->P0 == p1 && e->P1 == p0))
                return e;
        return nullptr;
    };
    Triangle *Network::findTriangle(Point *p0, Point *p1, Point *p2)
    {
        for (const auto &e : this->triangles)
            if ((e->P0 == p0 && e->P1 == p1 && e->P2 == p2) ||
                (e->P0 == p0 && e->P1 == p2 && e->P2 == p1) ||
                (e->P0 == p1 && e->P1 == p0 && e->P2 == p2) ||
                (e->P0 == p1 && e->P1 == p2 && e->P2 == p0) ||
                (e->P0 == p2 && e->P1 == p0 && e->P2 == p1) ||
                (e->P0 == p2 && e->P1 == p1 && e->P2 == p0))
                return e;
        return nullptr;
    };
    /* ------------------------------------------------------ */
}
Tddd ToX(geometric_network::Point *X) { return X->X; };

#endif