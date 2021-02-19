
#pragma once

#include "geom.h"
#include <algorithm>
#include <cstdint>
#include <string>
#include <cmath>

#include <map>
#include <unordered_map>
#include <set>
#include <list>
#include <vector>

#define logoutf1(x) printf(x "\n")
#define logoutf(format, ...) printf(format "\n", __VA_ARGS__)

struct PolyLine;
struct PolyDetector;

using CycleSet = std::set<uint32_t>;

enum class RmLinesType : int
{
    TakenTwice = 0, // rm lines taken twice
    Collinear, // dissolve collinear lines
    NoPointNeigh, // rm lines with one point not having neighbors
    PointConsumed, // point has 2 lines connected all taken
};

#if 1
using PointType = vec;
#else
// In case some additional members are needed
struct PointType : public vec
{
    PointType() {}
    
    PointType(vecType v) :
        vec(v)
    {}
    
    PointType(vecType xP, vecType yP, vecType zP):
        vec(xP, yP, zP)
    {}
    
    //bool ignore = false;
    uint32_t took = 0;
};
#endif

struct PolyCycle
{
    CycleSet idx;
    
    uint32_t startIdx, lastIdx;
    bool isClosed;
    bool fine;
    //uint32_t colStep = 0; // use it to increase the recursive colliniar length
    
    bool canBeClosed(PolyDetector &pd, uint32_t idToAdd) const;
    
    bool contains(uint32_t idP) const
    {
        return std::find(idx.begin(), idx.end(), idP) != idx.end();
    }
    
    uint32_t numCuts(PolyDetector &pd, const PolyLine &l) const;
    
    bool AddLineId(PolyDetector &pd, uint32_t id);
    
    std::string idxToString() const
    {
        std::string str;
        for (auto &n : idx)
        {
            str += std::to_string(n) + " ";
        }
        if (!str.empty())
            str.pop_back();
        return str;
    }
    
    std::string toString() const
    {
        std::string str = "C{[" + std::to_string(startIdx) + "] nLines:";
        str += std::to_string(idx.size());
        
        str += " [";
        str += idxToString();
        str += "]}";
        
        return str;
    }
    
    void print(const char *msg = nullptr) const
    {
        logoutf("%s%s", msg ? msg : "", toString().c_str());
    }
    
    bool Equals(const PolyCycle &p) const;
    
    bool convex(PolyDetector &pd) const;
    
    bool pointConsumed(PolyDetector &pd, uint32_t pid) const;
    
    bool accepted(PolyDetector &pd);
};
using PolyCycles = std::vector<PolyCycle>;

struct PolyLine
{
    PointType a, b;
    
    PolyLine() {}
    PolyLine(const PointType &aP, const PointType &bP):
        a(aP), b(bP)
    {}
    
    std::vector<uint32_t> intersections;
    std::set<uint32_t> intersectedLines;
    PointType center;
    uint32_t id = 0;
    uint32_t test0 = 0, test1 = 0;
    //obb::Line normal;
    int32_t origLine = -1;
    int32_t attr0 = -1; // used for testing
    
    bool ignore = false;
    uint32_t took = 0;
    uint32_t processed = 0;
    
    uint32_t aIdx = 0, bIdx = 0;
    
    uint32_t lastDissolveStep = 0;
    
    //obb::Line &calcNormal(PolyDetector &pd);
    
    void calcCenter();
    
    //bool HasCommonPoints(const PolyLine &line) const;
    bool HasCommonIdxPoints(const PolyLine &line) const;
    bool Equals(const PolyLine &line) const;
    
    void SortIntersectionsList(PolyDetector &pd);
    
    bool IntersectionPoint(const PolyLine &line, PointType &pos) const;
    bool LineLineIntersectionPoint(const PolyLine &line, PointType &pos) const;
    
    void CalculateFirstAndLastPoint();

    static bool HaveCommonPoints(const PolyLine &l1, const PolyLine &l2);
    static bool HaveCommonIdxPoints(const PolyLine &l1, const PolyLine &l2);
    
    static bool bCompareLineOrder(const PolyLine &l1, PolyLine &l2);
    static int iCompareLineOrder(const PolyLine &l1, PolyLine &l2);
    
    std::string neighToString(PolyDetector &pd, uint32_t *retNNeigh = nullptr) const;
    std::string toString(PolyDetector &pd) const;
    uint32_t numNeigh(PolyDetector &pd) const;
    uint32_t numIntersections(PolyDetector &pd) const;
    
    uint32_t canBeRemoved(PolyDetector &pd, RmLinesType type) const;
    
    PolyLine &mul(float m)
    {
        a.mul(m);
        b.mul(m);
        return *this;
    }
    PolyLine &add(const PointType &p)
    {
        a.add(p);
        b.add(p);
        return *this;
    }
    
    /*PolyLine &rot(float r, const vec axis = {0, 0, 1})
    {
        a.rotate(r, axis);
        b.rotate(r, axis);
        return *this;
    }*/
    
    void setIgnore(PolyDetector &pd, const char *msg);
    bool contains(const PointType &point) const;
    bool contains(const PolyLine &line) const;
    bool collinear(const PolyLine &l) const;
    //bool intersects(const PolyLine &line) const;
    
    uint32_t minPid() const
    {
        return std::min(aIdx, bIdx);
    }
    uint32_t maxPid() const
    {
        return std::max(aIdx, bIdx);
    }
    
    int32_t commonPid(const PolyLine &l) const
    {
        if (aIdx == l.aIdx || aIdx == l.bIdx) return aIdx;
        if (bIdx == l.aIdx || bIdx == l.bIdx) return bIdx;
        return -1;
    }
    
    uint32_t otherPid(uint32_t pid) const
    {
        return aIdx == pid ? bIdx : aIdx;
    }
    
    bool compareNeigh(PolyDetector &pd, uint32_t nid1, uint32_t nid2) const;
    bool sortNeigh(PolyDetector &pd) const;
    uint32_t &incTook(PolyDetector &pd);
    
    float angle(PolyDetector &pd, const PolyLine &l) const;
    bool betweenNeighbors(PolyDetector &pd, const PolyLine &l1, const PolyLine &l2) const;
};

struct PolyPol
{
    // members {
    std::vector<PointType> p;
    size_t firstIdx = 0;
    
    PolyCycle cycle;
    
    uint32_t id;
    uint32_t dissolveStep = 0;
    
    //obb::Polygon poly;
    double _area;
    // } members
    
    uint32_t GetCount() const { return (uint32_t)p.size(); }
    
    void CalculateFirstAndLastPoint();
    
    bool IsClosed() const;
    
    static bool HaveCommonVertex(const PolyPol &p1, const PolyPol &p2, size_t &i, size_t &j);
    bool Minus(const PolyPol &other);
    
    void addLine(const PolyLine &l);
    
    bool addPointChecked(const PointType &v);
    
    PointType center();
    
    void genTriangePoints(PolyDetector &pd);
    double TriangleArea(PolyDetector &pd);
    
    //obb::Polygon &updatePoly();
};

struct PolyDetector
{
    using LineVector = std::vector<PolyLine>;
    using LineList = std::list<PolyLine>;
    using PolyVector = std::vector<PolyPol>;
    
    // members {
    PolyCycles _cycles;
    
    //std::vector<CycleSet> processed;
    
    LineVector
        origLines, // from user
        lines; // active lines
    PolyVector polys;
    uint32_t verbose = 0;
    
    std::map<uint32_t, std::vector<uint32_t>> _neighbors; // key: lid, val: vec(neighborLine)
    std::map<uint32_t, std::vector<uint32_t>> collinearLineMap; // key: lid, val: vec(collinearLine)
    
    std::vector<PointType> intersectionPoints;
    std::map<uint32_t, std::vector<uint32_t>> pointToLines; // key: pid, val: list of nids
    //std::map<uint32_t, int32_t> lineIdToIdx; // key: lid, val: index in lines
    std::unordered_map<uint32_t, int32_t> lineIdToIdx; // key: lid, val: index in lines
    
    uint32_t dissolveCount = 0;
    
    // } members
    
    // ops {
    void reset();
    void AddLine(const PolyLine &line);
    bool DetectPolygons();
    PolyVector &getPolys() { return polys; }
    // } ops
    
    bool addPointToLine(uint32_t pid, uint32_t lid);
    
    PolyLine *findLine(uint32_t id, bool useIgnore = true);
    PolyLine *findLine(uint32_t pidA, uint32_t pidB, bool useIgnore = true);
    PolyLine *findOrigLine(uint32_t id);
    
    uint32_t GetPolyCount() const { return (uint32_t)polys.size(); };
    void SortLines();
    
    bool cycleProcessed(const PolyCycle &cycle) const;
    
    // The order of operations
    void RemoveZeroLengthLines();
    void RemoveOverlappings();
    void mergeOverlapped(PolyLine &lid1, PolyLine &lid2);
    uint32_t DetectAllIntersections();
    bool CreateLines();
    bool FindPolys();
    void SimplifyPolys(double smaller_polygon_length);
    
    PointType *getPointByID(size_t id);
    
    bool BuildCycle(uint32_t id, PolyCycle cycle); // as value!
    
    bool Overlap(const PolyCycle &c1, const PolyCycle &c2);
    
    PolyLine newLine(uint32_t i, uint32_t j, PolyLine &origLine);
    
    // collinearity
    void setCollinear(uint32_t l1, uint32_t l2);
    bool CollinearIdx(uint32_t l1, uint32_t l2);
    bool CollinearIdx(const PolyLine &l1, const PolyLine &l2);
    
    bool rmLines(RmLinesType type);
    bool rmEarPoints();
    
    // dissolve
    bool dissolveCollinearLine(PolyLine &l);
    bool dissolveCollinear(PolyLine &l1, PolyLine &l2);
    bool dissolve();
    
    void dumpLines(const char *msg, bool useIgnore = false);
    
    
};
