
#include "common.h"
#include "obb.h"
#include "PolyDetector.h"

#include <set>

static bool pointsDiffer(const PointType &a, const PointType &b)
{
    //return a.x != b.x || a.y != b.y;
    return fabs(a.x - b.x) > 1e-4f || fabs(a.y - b.y) > 1e-4f;
    //return a != b;
}

static int iComparePointOrder(const PointType &p1, const PointType &p2)
{
    if (p1.y < p2.y)
        return -1;
    else if (p1.y == p2.y)
    {
        if (p1.x < p2.x)
            return -1;
        else if (p1.x == p2.x)
            return 0;
    }
    // p1 is greater than p2
    return 1;
}

static bool bComparePointOrder(const PointType &p1, const PointType &p2)
{
    return iComparePointOrder(p1, p2) < 0;
}

bool PolyLine::bCompareLineOrder(const PolyLine &l1, PolyLine &l2)
{
    return iCompareLineOrder(l1, l2) < 0;
}
int PolyLine::iCompareLineOrder(const PolyLine &l1, PolyLine &l2)
{
    int result = iComparePointOrder(l1.a, l2.a);

    if (result == 0)
    {
        // in case lines share first point
        // we must order the lines by its slope

        auto dx1 = l1.b.x - l1.a.x;
        auto dy1 = l1.b.y - l1.a.y;
        auto dx2 = l2.b.x - l2.a.x;
        auto dy2 = l2.b.y - l2.a.y;

        // by definition of first and last point we are sure that dy > 0

        if (dx1>0 && dx2<0)
            // line 1 in 1st quadrant, line 2 in 2nd quadrant
            // this means line 2 cames first
            return 1;

        if (dx1<0 && dx2>0)
            // line 1 in 2nd quadrant, line 2 in 1st quadrant
            // this means line 1 cames first
            return -1;

        if (dx1 == 0) {
            // first line is vertical
            if (dx2>0)
                // second line in 1st quadrant
                // first line is previous
                return -1;

            if (dx2<0)
                // second line in 2nd quadrant
                // second line is previous
                return 1;
            // this should no happen
            return 0;
        }

        if (dx2 == 0) {
            // second line is vertical
            if (dx1>0)
                // first line in 1st quadrant
                // second line is previous
                return 1;

            if (dx1<0)
                // first line in 2nd quadrant
                // first line is previous
                return -1;

            // this should not happen
            return 0;
        }


        // calculate the slopes
        double m1 = dy1/dx1;
        double m2 = dy2/dx2;
        // line 1 and line 2 in 2nd quadrant
        if (m1 > m2)
            return -1;
        if (m1 < m2)
            return 1;
        
        // in this case we have the same slope in both lines,
        // which means that both lines are coincident.
        return 0;
    }

    return result;
}

bool RayCrossing(const PointType &rayOrigin, const PointType &p1, const PointType &p2)
{
    return (
        (p1.y > rayOrigin.y && p2.y <= rayOrigin.y)
        ||
        (p2.y > rayOrigin.y && p1.y <= rayOrigin.y)
    ) &&
    (p1.x+(rayOrigin.y-p1.y) / (p2.y-p1.y)*(p2.x-p1.x) < rayOrigin.x);
}

static float Area(const PointType &a, const PointType &b, const PointType &c)
{
    return 0.5 * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
}

static bool Collinear(const PointType &a, const PointType &b, const PointType &c)
{
    //return Area(a, b, c) == 0.0f;
    return fabs(Area(a, b, c)) <= 1e-3f;
}

//static bool Collinear(const PolyLine &l1, const PolyLine &l2)
//{
//    return (Collinear(l1.a, l2.a, l2.b) && Collinear(l1.b, l2.a, l2.b)) && l1.HasCommonPoints(l2);
//}

/***
* @desc indicates if this point is at left of the directed line from a to b
* @return true if is at left, false otherwise
* @see O'Rourke, Joseph, "Computational Geometry in C, 2nd Ed.", pp.29
*/
static bool Left(const PointType &a, const PointType &b, const PointType &c)
{
    return Area(a, b, c) > 0.0f;
}

/***
* @return true is this point is betwen a and b
* @note c must be collinear with a and b
* @see O'Rourke, Joseph, "Computational Geometry in C, 2nd Ed.", pp.32
*/
static bool Between(const PointType &p, const PointType &a, const PointType &b)
{
    // if this point is not collinear with a and b
    // then it cannot be between this two points
    if (! Collinear(p, a, b))
        return false;
    
    auto &_x = p.x;
    auto &_y = p.y;
    
    return
        ((a.x <= _x && _x <= b.x) && (a.y <= _y && _y <= b.y)) ||
        ((b.x <= _x && _x <= a.x) && (b.y <= _y && _y <= a.y));
}

/***
* @return true is this point is betwen a and b,
*        but is different from a and b
* @note c must be collinear with a and b
* @see O'Rourke, Joseph, "Computational Geometry in C, 2nd Ed.", pp.32
*/
static bool StrictBetween(const PointType &p, const PointType &a, const PointType &b)
{
    // first check if this point is between a and b
    if (!Between(p, a, b))
        return false;

    // if is between, lets check if it is not coincident with
    // one of them
    return pointsDiffer(a, p) && pointsDiffer(b, p);
}

static bool Overlap(const PolyLine &l1, const PolyLine &l2)
{
    auto &p1 = l1.a;
    auto &p2 = l1.b;
    auto &p3 = l2.a;
    auto &p4 = l2.b;
    
    // first see of all endpoints are collinear,
    // then see if any endpoint of one line lies on the other line
    return (Collinear(p1, p3, p4) && Collinear(p2, p3, p4)) &&
        ((l1.PolyContains(p3) || l1.PolyContains(p4)) ||
         (l2.PolyContains(p1) || l2.PolyContains(p2)));
}

/***
* @return a new simplified line if line_1 and line_2 overlaps, NULL otherwise
*/
int PolyLine::SimplifiedLine(const PolyLine &line_1, const PolyLine &line_2,
                                            PolyLine &ret)
{
    if (Overlap(line_1, line_2))
    {
        if (line_1.PolyContains(line_2))
        {
            ret = line_1;
            return 1;
        }
        if (line_2.PolyContains(line_1))
        {
            ret = line_2;
            return 2;
        }

        PointType new_line_start_point;
        PointType new_line_end_point;

        // detects which point of <line_1> must be removed
        if (Between(line_1.a, line_2.a, line_2.b)) {
            new_line_start_point = line_1.b;
        } else {
            new_line_start_point = line_1.a;
        }
        // detects which point of <line_2> must be removed
        if (Between(line_2.a, line_1.a, line_1.b)) {
            new_line_end_point = line_2.b;
        } else {
            new_line_end_point = line_2.a;
        }

        // create a new line
        ret = PolyLine(new_line_start_point, new_line_end_point);
        return 3;
    }

    return 0;

}

/***
* @return true if <line_1> contains <line_2>, false otherwise
*/
inline bool PolyLine::PolyContains(const PolyLine &line) const
{
    return PolyContains(line.a) && PolyContains(line.b);
}

/***
* @return true if <line_1> contains strictly <line_2>, false otherwise
*/
inline bool PolyLine::StrictContains(const PolyLine &line) const
{
    return StrictContains(line.a) && StrictContains(line.b);
}

inline bool PolyLine::PolyContains(const PointType &point) const
{
    return Between(point, a, b);
}


/***
* @return true if <line> contains strictly <point>, false otherwise
*/
inline bool PolyLine::StrictContains(const PointType &point) const
{
    return StrictBetween(point, a, b);
}

bool PolyLine::HasCommonPoints(const PolyLine &line) const
{
    return HaveCommonPoints(*this, line);
}

bool PolyLine::HaveCommonPoints(const PolyLine &l1, const PolyLine &l2)
{
    return (!pointsDiffer(l1.a, l2.a) || !pointsDiffer(l1.b, l2.b) ||
            !pointsDiffer(l1.a, l2.b) || !pointsDiffer(l1.b, l2.a));
}
/**
* param a pointer to a line
* @return true if line intersects with this, false otherwise
*/
bool PolyLine::PolyIntersects(const PolyLine &line) const
{
    if (PolyIntersectsProper(line))
        return true;
    
    auto &c = line.a;
    auto &d = line.b;
    
    return (Between(c, a, b) || Between(d, a, b) || Between(a, c, d) || Between(b, c, d));
}

/***
* param a pointer to a line
* @return true is line have a proper intersection with this, false otherwise
* @note a proper intersection occurs when two segments intersects at a point
*      interior to both
* @see O'Rourke, Joseph, "Computational Geometry in C, 2nd Ed.", pp.30
*/
bool PolyLine::PolyIntersectsProper(const PolyLine& line) const
{
    auto &c = line.a;
    auto &d = line.b;
    
    // Eliminate improper cases
    if (Collinear(c, a, b) || Collinear(d, a, b) || Collinear(a, c, d) || Collinear(b, c, d))
        return false;

    return (Left(c, a, b) ^ Left(d, a, b)) && (Left(a, c, d) ^ Left(b, c, d) );
}

bool PolyLine::IntersectionPoint(const PolyLine &line, PointType &pos) const
{
    // in case two lines share just vertices it was no intersection
    if (!pointsDiffer(line.a, a) ||
        !pointsDiffer(line.a, b) ||
        !pointsDiffer(line.b, a) ||
        !pointsDiffer(line.b, b))
    {
        return false;
    }
    
    double s, num, denom;
    
    auto &c = line.a;
    auto &d = line.b;

    // first we calculate the denominator of equation
    denom = a.x * (d.y - c.y) +
        b.x * (c.y - d.y) +
        d.x * (b.y - a.y) +
        c.x * (a.y - b.y);

    // see if segments are parallel
    if (denom == 0.0)
    {
        logoutf("zero denom!");
        return false;
    }

    num = a.x * (d.y - c.y) +
        c.x * (a.y - d.y) +
        d.x * (c.y - a.y);

    s = num / denom;
/*
    num = c.x * ( c.y - b.y) +
        b.x * ( a.y - c.y) +
        c.x * ( b.y - a.y);

    t = num / denom;
*/
    
    pos = PointType(a.x + s * (b.x - a.x), a.y + s * (b.y - a.y), 0.0f);
    
    return true;
}

void PolyLine::SortIntersectionsList()
{
    std::sort(intersections.begin(), intersections.end(), bComparePointOrder);
}
/***
* @desc Removes intersections between lines
* @return true if successfull, false otherwise
* @see Bentley-Ottmann Algorithm
*/
bool PolyDetector::RemoveIntersections()
{
    logoutf("Line intersection removal");
    
    // prior to removing overlapping, one must
    // remove all zero length line, otherwise the results
    // will be unpredictable
    RemoveZeroLengthLines();

    // then we must remove line overlapping in order to run
    // the Bentley-Ottmann Algorithm
    RemoveOverlappings();

    // finally we detect intersections between lines
    //int intersection_count = DetectIntersections();
    int intersection_count = DetectAllIntersections();
    if (!silent)
        logoutf("Detected %d intersections", intersection_count);
    
#if 0
    for (auto &line : lines)
        logoutf("line nIntersections:%u", uint32_t(line.intersections.size()));
    exit(0);
#endif

    LineVector createdLines;
    //ITERATE_ALL_ENTITIES(_lines_array,ResetFlag(FLAG_REMOVE));
    for (auto &line: lines)
        line.toRemove = false;

    // sweep all lines
    for (auto &line : lines)
    {
        // check if current line has intersections
        if (!line.intersections.empty())
        {
            auto point = line.a;
            line.SortIntersectionsList();
            //logoutf("line nIntersections:%u", uint32_t(line.intersections.size()));
            for (auto &intersection : line.intersections)
            {
                createdLines.push_back(PolyLine(point, intersection));
                point = intersection;
            }
            createdLines.push_back(PolyLine(point, line.b));
            
            line.toRemove = true;
        }
    }

    // Now lets remove the flagged lines
    for (auto it = lines.begin(); it != lines.end(); )
    {
        if (it->toRemove)
            it = lines.erase(it);
        else
            ++it;
    }
    for (auto &l : createdLines)
    {
        lines.push_back(l);
    }
    
    logoutf("%u lines after intersection removal", uint32_t(lines.size()));
    
    return true;
}

/***
* @desc removes all lines with zero length
*/
void PolyDetector::RemoveZeroLengthLines(void)
{
    for (auto it = lines.begin(); it != lines.end(); )
    {
        if (!pointsDiffer(it->a, it->b)) // find a zero length line
            it = lines.erase(it);
        else
            ++it;
    }
}

/***
* @descr removes line overlappings
* @note must be called before applying Bentley-Ottmann algorithm
*/
void PolyDetector::RemoveOverlappings()
{
    size_t i, j, count = lines.size();
    
    PolyLine line;

    // lets find overlapping lines
    for (i = 0; i < count; i++)
    {
        auto &line_i = lines[i];

        for(j = i + 1; j < count; j++)
        {
            auto &line_j = lines[j];
            if (::Overlap(line_i, line_j))
            {
                int ret = PolyLine::SimplifiedLine(line_i, line_j, line);

                if (ret == 1)
                {
                    // must remove line_j
                    lines.erase(lines.begin() + j);
                    j--;
                    count--;
                }
                else
                {
                    if (ret != 2)
                    {
                        // must remove both line_i and line_j and add a new one
                        lines.erase(lines.begin() + j);
                        lines.push_back(line);
                    }

                    // must remove line_i
                    lines.erase(lines.begin() + i);

                    // update counters
                    i--;
                    count --;

                    // skip inner loop an go to next step of outer loop
                    break;
                }
            }

        }
    }
}

uint32_t PolyDetector::DetectAllIntersections()
{
    uint32_t ret = 0;
    size_t counter = lines.size();
    PointType intersection;
    for (size_t i = 0; i < counter; i++)
    {
        auto &l1 = lines[i];
        for (size_t j = i + 1; j < counter; j++)
        {
            auto &l2 = lines[j];
            // if the current line intersects with other line
            if (l1.PolyIntersects(l2))
            {
                if (l1.IntersectionPoint(l2, intersection))
                {
                    l1.intersections.push_back(intersection);
                    l2.intersections.push_back(intersection);
                    ret++;
                }
            }
        }
    }
    return ret;
}

void PolyLine::CalculateFirstAndLastPoint()
{
    if (!bComparePointOrder(a, b))
    {
        std::swap(a.x, b.x);
        std::swap(a.y, b.y);
        std::swap(a.z, b.z);
    }
}

/***
* @descr sort the lines
*/
void PolyDetector::SortLines(void)
{
    for (auto &line : lines)
        line.CalculateFirstAndLastPoint();
    std::sort(lines.begin(), lines.end(), PolyLine::bCompareLineOrder);
}

/**
* @descr Performs polygon detection
*/
bool PolyDetector::DetectPolygons()
{
    logoutf("Polygon detection");

    if (!silent)
        logoutf("Line set contains %d lines.", GetLineCount());

    if (!RemoveIntersections())
    {
        logoutf("Error: Could not successfully remove line intersections.");
        return false;
    }

    if (!silent)
        logoutf("After removal, line set contains %d lines.", GetLineCount());

    SortLines();
    
    if (!FindPolys())
    {
        logoutf("Error constructing the polygon set.");
        return false;
    }

    //if (!silent)
        logoutf("Polygon set contains %d polygons.", GetPolyCount());
    
    //SimplifyPolys(0.0);
    //logoutf("Polygon set contains %d polygons after simplification!.", GetPolyCount());

    return true;
}

bool PolyPol::IsClosed() const
{
    return !p.empty() && !pointsDiffer(p[0], p.back());
}

PointType PolyPol::center()
{
    if (p.empty())
        return PointType(0);
    
    PointType c(0);
    for (auto &pt : p)
        c.add(pt);
    c.div(p.size());
    
    return c;
}

/***
* @return true if polylines have a common vertex, false otherwise
* param p1, p2 pointers to polylines
* param i, j pointers to the indices of first common vertices found
*/
bool PolyPol::HaveCommonVertex(const PolyPol &p1, const PolyPol &p2, size_t &i, size_t &j)
{
    // sweeps all vertices in p1
    for (i = 0; i < p1.p.size(); ++i)
    {
        auto &v1 = p1.p[i];
        
        for (j = 0; j < p2.p.size(); ++j)
        {
            auto &v2 = p2.p[j];
            if (!pointsDiffer(v1, v2))
                return true;
        }
    }

    return false;
}

/***
* @desc simplifies current polygon, subtracting <p>, in case there are
*       an single relationship between them
* @return true if polygon were changed, false otherwise
*/
bool PolyPol::Minus(const PolyPol &other)
{
    for (auto itOther = other.p.begin(); itOther != other.p.end(); ++itOther)
    {
        for (auto it = p.begin(); it != p.end();)
        {
            if (!pointsDiffer(*itOther, *it))
            {
                it = p.erase(it);
            }
            else
                ++it;
        }
    }
    return true;
}

void PolyPol::CalculateFirstAndLastPoint()
{
    // if there are only one vertex it is not a polyline
    if (p.size() < 2)
    {
        p.clear();
        return;
    }
    
    if (!IsClosed())
    {
        p.clear();
        return;
    }
    
    // the case of the closed polyline
    // here we're going to find the first point by seeing them all
    //Point2D *vertex=NULL, *first_vertex=_vertex_array[0];
    PointType *firstVertex = nullptr;
    for (size_t i = 0; i < p.size(); ++i)
    {
        auto &vertex = p[i];
        if (!firstVertex || !bComparePointOrder(*firstVertex, vertex))
        {
            firstVertex = &vertex;
            firstIdx = i;
        }
    }
}

bool PolyPol::IsAdjacent(const PolyPol &other, bool strict) const
{
    PointType previous_v1 = p.back();
    PointType previous_v2 = other.p.back();

    for (auto &v1 : p)
    {
        PolyLine l1(previous_v1, v1);
        
        for (auto &v2 : other.p)
        {
            PolyLine l2(previous_v2, v2);

            //if (strict ? StrictOverlap(l1, l2) : Overlap(l1, l2))
            if (Overlap(l1, l2))
            {
                return true;
            }
            
            previous_v2 = v2;
        }

        previous_v1 = v1;
    }
    
    return false;
}

bool PolyPol::PolyContains(const PointType &point, bool strict)
{
    // in case of an open polyline the point cannot be inside
    if (!IsClosed()) return false;

    bool inside = false;
    size_t vertex_count = p.size();

    if (vertex_count > 0)
    {
        PointType *previous_vertex = nullptr;
        
        for (size_t i = 0; i < vertex_count - 1; ++i)
        {
            auto &current_vertex = p[i];
            
            if (previous_vertex)
            {
                // let's see if point is lays on the edge
                if (Between(point, current_vertex, *previous_vertex))
                    return true;

                if (RayCrossing(point, *previous_vertex, current_vertex))
                    inside = !inside;
            }

            previous_vertex = &current_vertex;
        }
    }

    return inside;
}
        
bool PolyPol::PolyContains(const PolyPol &other, bool strict)
{
    // first lets see if all vertices in polyline lay inside or in the border of the polygon
    for (auto &v : other.p)
    {
        if (!obb::Polygon::Contains(v))
            return false;
    }
    
    return true;
}

/***
* @desc simplifies the polygons in this set
* @note removes inclusions and disposes small polygons
*/
void PolyDetector::SimplifyPolys(double smaller_polygon_length)
{
    logoutf("Polygon set simplification");

    // remove small polygons
    uint32_t nRemoved = 0;
    for (auto it = polys.begin(); it != polys.end(); )
    {
        if (it->GetCount() < smaller_polygon_length)
        {
            it = polys.erase(it);
            nRemoved++;
        }
        else
            ++it;
    }

    if (nRemoved)
        logoutf("Removed %d small polygons", nRemoved);
    
#if 0
    // check adjacencies
    nRemoved = 0;
    for (size_t i = 0; i < polys.size(); i++)
    {
        auto &p1 = polys[i];
        for (size_t j = i + 1; j < polys.size(); j++)
        {
            auto &p2 = polys[j];

            // see if p1 and p2 are strictly adjacent
            if (p1.IsAdjacent(p2))
            {
                // check if p1 contains p2
                //if (p1.Polygon::Contains(p2))
                if (p1.PolyContains(p2))
                {
                    p1.Minus(p2);
                    //p2.p.clear();
                    nRemoved++;
                }
                //else if (p2.Polygon::Contains(p1))
                else if (p2.PolyContains(p1))
                {
                    // just when p1 does not contains p2
                    // we need to check if p2 contains p1
                    p2.Minus(p1);
                    //p1.p.clear();
                    nRemoved++;
                }
            }
        }
    }
    //if (!silent)
    if (nRemoved)
        logoutf("Merged %d contained polygons", nRemoved);
#endif
}

void PolyDetector::AddLine(const obb::LineSegment &line)
{
    AddLine(PolyLine(line.a, line.b));
}

void PolyDetector::AddLine(const PolyLine &line)
{
    lines.push_back(line);
}

void PolyPol::addLine(const PolyLine &l)
{
    p.push_back(l.a);
    p.push_back(l.b);
}

PolyLine *PolyDetector::findLine(uint32_t id)
{
    for (auto &l : lines)
    {
        if (l.id == id)
            return &l;
    }
    logoutf("Cannot find line with id %u. nlines:%u", id, uint32_t(lines.size()));
    return nullptr;
}

bool similarCycle(const PolyCycles &cycles, const PolyCycle &cycle)
{
    for (auto &c : cycles)
    {
        bool equal = true;

        // set already sorted
        for (auto itC = c.idx.cbegin(), itCycle = cycle.idx.cbegin(); itC != c.idx.cend();)
        {
            if (*itC != *itCycle)
            {
                equal = false;
                break;
            }
            bool cEnd = ++itC == c.idx.cend();
            bool cycleEnd = ++itCycle == cycle.idx.cend();
            if (cEnd || cycleEnd)
            {
                if (cEnd != cycleEnd) // contained!
                {
                    return true;
                }
                break;
            }
        }
        
        if (equal)
            return true;
    }
    return false;
}

bool PolyDetector::Overlap(const PolyCycle &c1, const PolyCycle &c2)
{
    for (auto &n1 : c1.idx)
    {
        for (auto &n2 : c2.idx)
        {
            if (n1 != n2)
            {
                auto l1 = findLine(n1);
                auto l2 = findLine(n1);
                if (l1 && l2)
                {
                    if (::Overlap(*l1, *l2))
                        return true;
                }
            }
        }
    }
    return false;
}

bool PolyCycle::AddLineId(PolyDetector &pd, uint32_t id)
{
    auto l = pd.findLine(id);
    if (!l) return false;
    
    l->test0 = 0; // a cnt
    l->test1 = 0; // b cnt
    
    for (auto &id1 : idx)
    {
        auto l1 = pd.findLine(id1);
        if (!l1) continue;
        
        if (!pointsDiffer(l->a, l1->a) || !pointsDiffer(l->a, l1->b))
        {
            if (Collinear(l->a, l1->a, l1->b) && Collinear(l->b, l1->a, l1->b)) // collinear points
                return false;
            l->test0++;
        }
        if (!pointsDiffer(l->b, l1->a) || !pointsDiffer(l->b, l1->b))
        {
            if (Collinear(l->a, l1->a, l1->b) && Collinear(l->b, l1->a, l1->b)) // collinear points
                return false;
            l->test1++;
        }
        
        if (l->test0 >= 2 || l->test1 >= 2)
        {
            //logoutf("line %u can't be added to cycle! ta:%d tb:%d ids: [%s]", id, l->test0, l->test1, idxToString().c_str());
            return false;
        }
        
        obb::LineSegment ls(l->center, l1->center);
        for (auto &id2 : idx)
        {
            if (id2 != id && id2 != id1)
            {
                auto l2 = pd.findLine(id2);
                if (l2 && ls.Intersects(*l2))
                {
                    return false;
                }
            }
        }
    }
    idx.insert(id);
    return true;
}

static bool equalCycles(const CycleSet &a, const CycleSet &b)
{
    if (a.size() != b.size())
        return false;
    for (auto aIt = a.begin(), bIt = b.begin(); aIt != a.end(); ++aIt, ++bIt)
        if (*aIt != *bIt)
            return false;
    return true;
}

bool PolyDetector::CycleProcessed(const CycleSet &cycle) const
{
    for (auto &c : processed)
        if (equalCycles(c, cycle))
            return true;
    return false;
}

bool PolyDetector::BuildCycle(uint32_t id, PolyCycle cycle) // as value!
{
    //if (cycle.startIdx != id && cycle.idx.size() > 2)
    //    return true;
    
    //if (cycle.idx.size() >= 15)
    //    return true;
    
    auto l = findLine(id);
    if (!l)
        return true;
    
    if (_neighbors[id].size() < 4)
    {
        return true;
    }
    
    if (!silent)
        cycle.print((std::string("[PROC:") + std::to_string(id) + std::string("]")).c_str());
    
    if (cycle.canBeClosed(id))
    {
        cycle.isClosed = true;
        
        if (!silent)
            cycle.print("[CLOSED] ");
        
        if (!similarCycle(_cycles, cycle))
        {
            _cycles.push_back(cycle);
            if (!silent)
                cycle.print("[ACCEPTED] ");
            //logoutf("nCycles:%u", uint32_t(_cycles.size()));
        }
        else
        {
            if (!silent)
                cycle.print("[CYCLE_EXISTS] ");
        }
        
        return true;
    }
    
    if (!cycle.AddLineId(*this, l->id))
    {
        return true;
    }
    
    if (CycleProcessed(cycle.idx))
    {
        if (!silent)
            cycle.print("[PROCESSED!] ");
        return true;
    }
    processed.insert(cycle.idx);
    
    //if (polyIndices primaryId == polyIndices[0])

    for (auto &nid : _neighbors[id])
    {
        if (_neighbors[nid].size() < 4) continue;
        if (cycle.canBeClosed(nid) || !cycle.contains(nid))
        {
            if (!silent)
                cycle.print((std::string("[NEIGN:") + std::to_string(nid) + std::string("]")).c_str());
            BuildCycle(nid, cycle);
        }
    }
    
    return true;
}

bool PolyDetector::FindPolys()
{
    logoutf("Polygon set construction (MY)");
    
    //silent = 0;
    
    processed.clear();
    polys.clear();
    
    uint32_t n = 0;
    // assign line ids
    for (auto &l : lines)
    {
        l.id = n++;
    }
    
    // Build neighbors
    std::vector<uint32_t> neigh;
    _neighbors.clear();
    for (uint32_t i = 0; i < lines.size(); ++i)
    {
        auto &l1 = lines[i];
        neigh.clear();
        for (uint32_t j = i + 1; j < lines.size(); ++j)
        {
            auto &l2 = lines[j];
            if (l1.HasCommonPoints(l2))
            {
                neigh.push_back(j);
            }
        }
        if (!neigh.empty())
        {
            std::sort(neigh.begin(), neigh.end(), [this, &i, &l1](const uint32_t &a, const uint32_t &b) {
                auto da = lines[a].center.squaredist(l1.center);
                auto db = lines[b].center.squaredist(l1.center);
                return da < db;
            });
            for (auto &n : neigh)
            {
                _neighbors[i].push_back(n);
                _neighbors[n].push_back(i);
            }
        }
    }
    logoutf("neighbors MatSize:%u", uint32_t(_neighbors.size()));
    
    if (!silent)
    {
        for (auto &kv : _neighbors)
        {
            std::string str;
            for (auto &n : _neighbors[kv.first])
            {
                str += std::to_string(n) + " ";
            }
            logoutf("[%u] nNeigh:%u [%s]", kv.first, uint32_t(_neighbors[kv.first].size()), str.c_str());
        }
    }
    
    uint32_t step = 0;
    for (auto &kv : _neighbors) // point by point
    {
        step++;
        
        if (_neighbors[kv.first].size() < 4) continue;
        
        if (!silent)
            logoutf("----------- BEGIN %u {", kv.first);
        PolyCycle cycle;
        cycle.startIdx = kv.first;
        BuildCycle(kv.first, cycle);
        if (!silent)
            logoutf("} ----------- END %u", kv.first);
    }
    
    std::sort(_cycles.begin(), _cycles.end(), [](const PolyCycle &a, const PolyCycle &b) {
        return a.idxToString().compare(b.idxToString()) < 0;
    });
    
#if 0
    for (uint32_t i = 0; i < _cycles.size(); ++i)
    {
        auto &c1 = _cycles[i];
        if (!c1.fine) continue;
        for (uint32_t j = 0; j < _cycles.size(); ++j)
        {
            if (i == j) continue;
            auto &c2 = _cycles[j];
            if (!c2.fine) continue;
            
            if (/* DISABLES CODE */ (0))
            {
                if (c2.idx.size() > c1.idx.size())
                {
                    bool included = true;
                    for (auto &n1 : c1.idx)
                    {
                        if (std::find(c2.idx.begin(), c2.idx.end(), n1) == c2.idx.end())
                        {
                            included = false;
                            break;
                        }
                    }
                    if (included)
                    {
                        c2.fine = false;
                        //c1.fine = false;
                    }
                }
            }
        }
    }
#endif
    
    for (auto &cycle : _cycles)
    {
        //if (!cycle.fine) continue;
        
        if (!silent)
            cycle.print(("[FINI F:" + std::to_string(int(cycle.fine)) + "] ").c_str());
        
        if (!cycle.fine) continue;
        
        PolyPol poly;
        PointType last(0);
        for (auto &id : cycle.idx)
        {
            auto lPtr = findLine(id);
            if (lPtr)
            {
                lPtr->CalculateFirstAndLastPoint();
                
                //poly.addLine(*lPtr);
                poly.p.push_back(lPtr->a);
                
                last = lPtr->b;
            }
        }
        if (!cycle.idx.empty())
            poly.p.push_back(last);
        
#if 1
        uint32_t nRemoved = 0;
        for (int i = 0; i < poly.p.size(); ++i)
        {
            for (uint32_t j = i + 1; j < poly.p.size(); ++j)
            {
                if (!pointsDiffer(poly.p[i], poly.p[j]))
                {
                    poly.p.erase(poly.p.begin() + i);
                    i--;
                    nRemoved++;
                    break;
                }
            }
        }
        if (nRemoved)
        {
            //logoutf("Removed %u poly pts. poly pts size:%u", nRemoved, uint32_t(poly.p.size()));
        }
#endif
        poly.id = cycle.startIdx;
        poly.cycle = cycle;
        
        polys.push_back(std::move(poly));
    }
    
    return true;
}

obb::Line &PolyLine::calcNormal(PolyDetector &pd)
{
    auto line = pd.findLine(id);
    if (!line)
    {
        return normal;
    }
    vec n = vec(line->a).sub(line->b).cross(vec(0,0,1));
    n.normalize();
    normal = obb::Line(line->center, n);
    return normal;
}

void PolyPol::genTriangePoints(PolyDetector &pd)
{
    std::vector<size_t> ret;
    
    std::vector<PolyLine *> lines;
    for (auto &n : cycle.idx)
    {
        auto l = pd.findLine(n);
        if (!l) continue;
        lines.push_back(l);
    }
    
    for (auto &l : lines)
    {
        l->calcNormal(pd);
        for (auto &l1 : lines)
        {
            if (l != l1)
            {
                if (l->Intersects(*l1))
                {
                    vec pt = l->ClosestPoint(*l1);
                    pt.y -= .01f;
                    p.push_back(pt);
                    break;
                }
            }
        }
    }
}

double PolyPol::TriangleArea(PolyDetector &pd)
{
    if (p.size() <= 2)
        return 0.0f;
    double area = 0.0;

    auto c = center();
    for (auto &n : cycle.idx)
    {
        auto l = pd.findLine(n);
        if (!l)
            return 0.0f;
        area += abs(::Area(c, l->a, l->b));
    }
    //logoutf("area:%f", area);
    return area;
}
