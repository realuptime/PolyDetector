/*
MIT License

Copyright (c) 2021 Catalin Stan (catalinux)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#pragma once

#include "geom.h"

#include <algorithm>
#include <cstdint>
#include <string>
#include <map>
#include <set>
#include <vector>

#define logoutf(format, ...) printf(format "\n", __VA_ARGS__)

struct PolyDetector;

using PointType = vec; // define whatever you want in here
using CycleSet = std::set<uint32_t>;

struct PolyCycle
{
    CycleSet idx;
    
    uint32_t startIdx = 0;
    bool isClosed = false;
    bool fine = true;
    
    bool canBeClosed(uint32_t idToAdd) const
    {
        return idx.size() > 2 && startIdx == idToAdd;
    }
    
    bool contains(uint32_t idP) const
    {
        return std::find(idx.begin(), idx.end(), idP) != idx.end();
    }
    
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
        std::string str = "[" + std::to_string(startIdx) + "] nLines:";
        str += std::to_string(idx.size());
        
        str += " [";
        str += idxToString();
        str += "]";
        
        return str;
    }
    
    void print(const char *msg = nullptr) const
    {
        logoutf("%s%s", msg ? msg : "", toString().c_str());
    }
};
using PolyCycles = std::vector<PolyCycle>;

struct PolyLine 
{
    vec a, b;

    PolyLine() {}
    PolyLine(const PointType &aP, const PointType &bP): a(aP), b(bP)
    {
        CalculateFirstAndLastPoint();
        center = a;
        center.add(b);
        center.mul(0.5f);
    }
    
    bool toRemove = false;
    std::vector<PointType> intersections;
    PointType center;
    uint32_t id = 0;
    uint32_t test0 = 0, test1 = 0;

    //obb::Line normal;
    //obb::Line &calcNormal(PolyDetector &pd);
    
    bool HasCommonPoints(const PolyLine &line) const;
    
    void SortIntersectionsList();
    inline bool HaveIntersections() { return !intersections.empty(); };
    
    inline bool PolyContains(const PolyLine &line) const;
    inline bool PolyContains(const PointType &point) const;
    
    inline bool StrictContains(const PolyLine &line) const;
    inline bool StrictContains(const PointType &line) const;
    
    inline bool PolyIntersects(const PolyLine &line) const;
    inline bool PolyIntersectsProper(const PolyLine &line) const;
    
    bool IntersectionPoint(const PolyLine &line, PointType &pos) const;
    
    void CalculateFirstAndLastPoint();

    static int SimplifiedLine(const PolyLine &l1, const PolyLine &l2, PolyLine &ret);
    static bool HaveCommonPoints(const PolyLine &l1, const PolyLine &l2);
    
    static bool bCompareLineOrder(const PolyLine &l1, PolyLine &l2);
    static int iCompareLineOrder(const PolyLine &l1, PolyLine &l2);
};

struct PolyPol
{
    std::vector<vec> p;

    uint32_t GetCount() const { return (uint32_t)p.size(); };
    
    bool IsAdjacent(const PolyPol &p, bool strict = false) const;
    
    bool PolyContains(const PolyPol &other, bool strict = false);
    bool PolyContains(const PointType &point, bool strict = false);
    
    void CalculateFirstAndLastPoint();
    
    bool IsClosed() const;
    
    static bool HaveCommonVertex(const PolyPol &p1, const PolyPol &p2, size_t &i, size_t &j);
    bool Minus(const PolyPol &other);
    
    void addLine(const PolyLine &l);
    
    PointType center();
    
    void genTriangePoints(PolyDetector &pd);
    double TriangleArea(PolyDetector &pd);
    double _area;
    
    size_t firstIdx = 0;
    
    PolyCycle cycle;
    
    uint32_t id;
};

struct PolyDetector
{
    std::set<CycleSet> processed;
    using LineVector = std::vector<PolyLine>;
    using PolyVector = std::vector<PolyPol>;

    LineVector lines;
    PolyVector polys;
    bool silent = true;
    
    uint32_t GetLineCount() const { return (uint32_t)lines.size(); };
    uint32_t GetPolyCount() const { return (uint32_t)polys.size(); };
    void SortLines();
    
    void AddLine(const PolyLine &line);
    
    bool CycleProcessed(const CycleSet &cycle) const;
    
    bool RemoveIntersections();
    
    void RemoveZeroLengthLines();
    void RemoveOverlappings();
    uint32_t DetectAllIntersections();
    bool DetectPolygons();
    
    PointType *getPointByID(size_t id);
    
    bool BuildCycle(uint32_t id, PolyCycle cycle); // as value!
    bool FindPolys();
    
    // Simplifies the polygon set
    void SimplifyPolys(double smaller_polygon_length);
    
    PolyCycles _cycles;
    std::map<uint32_t, std::vector<uint32_t>> _neighbors;
    PolyLine *findLine(uint32_t id);
    
    bool Overlap(const PolyCycle &c1, const PolyCycle &c2);
};
