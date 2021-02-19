
#include <assert.h>
#include "PolyDetector.h"

#include <set>

#if 0
#undef logoutf
#define logoutf(...) do {} while(false)
#endif

#define arToStr(arg) #arg

//static const float minPointDiff = .1f;
//static const float minPointDiff = 1e-2f;
//static const float minPointDiff = 1e-3f;
//static const float minPointDiff = 1e-4f;
static const float minPointDiff = 1e-5f;

static const float minPointDiffSq = minPointDiff * minPointDiff;

const char *RmLinesTypeStr(RmLinesType type)
{
    switch (type)
    {
        case RmLinesType::TakenTwice: return "TakenTwice";
        case RmLinesType::Collinear: return "Collinear";
        case RmLinesType::NoPointNeigh: return "NoPointNeigh";
        case RmLinesType::PointConsumed: return "PointConsumed";
    }
    return "UNKN";
}

// Given three colinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
static bool onSegment(const PointType &p, const PointType &q, const PointType &r)
{
    if (q.x <= std::max(p.x, r.x) && q.x >= std::min(p.x, r.x) &&
        q.y <= std::max(p.y, r.y) && q.y >= std::min(p.y, r.y))
       return true;
  
    return false;
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
static int orientation(const PointType &p, const PointType &q, const PointType &r)
{
    // See https://www.geeksforgeeks.org/orientation-3-ordered-points/
    // for details of below formula.
    double val =
            (q.y - p.y) * (r.x - q.x) -
            (q.x - p.x) * (r.y - q.y);
  
    if (val == 0.0) return 0;  // colinear
    //if (abs(val) <= 1e-7f) return 0;  // colinear
  
    return (val > 0) ? 1: 2; // clock or counterclock wise
}

static bool collinearVecs(const PointType &p, const PointType &q, const PointType &r)
{
    return orientation(p, q, r) == 0;
}
  
// The main function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
// https://www.cdn.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
static bool doIntersect(const PointType &p1, const PointType &q1, const PointType &p2, const PointType &q2)
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);
  
    // General case
    if (o1 != o2 && o3 != o4)
        return true;
  
    // Special Cases
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;
  
    // p1, q1 and q2 are colinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;
  
    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;
  
     // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;
  
    return false; // Doesn't fall in any of the above cases
}

static bool pointsDiffer(const PointType &a, const PointType &b, bool aprox = true)
{
    // max precision is mandatory since this can break convex polys!
    if (aprox)
        return a.squaredist(b) >= minPointDiffSq;
    return a.x != b.x || a.y != b.y;
}

/***
* @return true is this point is betwen a and b
* @note c must be collinear with a and b
* @see O'Rourke, Joseph, "Computational Geometry in C, 2nd Ed.", pp.32
*/
static bool between(const PointType &p, const PointType &a, const PointType &b)
{
    // if this point is not collinear with a and b
    // then it cannot be between this two points
    if (!collinearVecs(p, a, b))
        return false;
    
    auto &_x = p.x;
    auto &_y = p.y;
    
    return
        ((a.x <= _x && _x <= b.x) && (a.y <= _y && _y <= b.y)) ||
        ((b.x <= _x && _x <= a.x) && (b.y <= _y && _y <= a.y));
}

bool PolyLine::contains(const PolyLine &line) const
{
    return contains(line.a) && contains(line.b);
}

bool PolyLine::contains(const PointType &point) const
{
    return between(point, a, b);
}

bool PolyLine::collinear(const PolyLine &line) const
{
    return !doIntersect(a, b, line.a, line.b);
}

//bool PolyLine::intersects(const PolyLine &line) const
//{
//    return doIntersect(a, b, line.a, line.b);
//}

bool PolyLine::IntersectionPoint(const PolyLine &line, PointType &pos) const
{
    return LineLineIntersectionPoint(line, pos);
}

static bool overlap(const PolyLine &l1, const PolyLine &l2)
{
    return (collinearVecs(l1.a, l2.a, l2.b) && collinearVecs(l1.b, l2.a, l2.b)) &&
        ((l1.contains(l2.a) || l1.contains(l2.b)) ||
         (l2.contains(l1.a) || l2.contains(l1.b)));
}

/***
* @return a new simplified line if line_1 and line_2 overlaps, NULL otherwise
*/
static int simplifiedLine(const PolyLine &line_1, const PolyLine &line_2, PolyLine &ret)
{
    if (overlap(line_1, line_2))
    {
        if (line_1.contains(line_2))
        {
            ret = line_1;
            return 1;
        }
        if (line_2.contains(line_1))
        {
            ret = line_2;
            return 2;
        }

        PointType new_line_start_point;
        PointType new_line_end_point;

        // detects which point of <line_1> must be removed
        if (between(line_1.a, line_2.a, line_2.b)) {
            new_line_start_point = line_1.b;
        } else {
            new_line_start_point = line_1.a;
        }
        // detects which point of <line_2> must be removed
        if (between(line_2.a, line_1.a, line_1.b)) {
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

static float Area(const PointType &a, const PointType &b, const PointType &c)
{
    return 0.5 * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
}

/*
bool PolyLine::Equals(const PolyLine &line) const
{
    return
        (!pointsDiffer(line.a, a) || !pointsDiffer(line.a, b)) &&
        (!pointsDiffer(line.b, a) || !pointsDiffer(line.b, b));
}
*/

bool PolyLine::HasCommonIdxPoints(const PolyLine &line) const
{
    return
        aIdx == line.aIdx ||
        aIdx == line.bIdx ||
        bIdx == line.aIdx ||
        bIdx == line.bIdx;
}

/*
bool PolyLine::HasCommonPoints(const PolyLine &line) const
{
    return HaveCommonPoints(*this, line);
}

bool PolyLine::HaveCommonPoints(const PolyLine &l1, const PolyLine &l2)
{
    return (!pointsDiffer(l1.a, l2.a) || !pointsDiffer(l1.b, l2.b) ||
            !pointsDiffer(l1.a, l2.b) || !pointsDiffer(l1.b, l2.a));
}
*/

// https://www.geeksforgeeks.org/program-for-point-of-intersection-of-two-lines/
bool PolyLine::LineLineIntersectionPoint(const PolyLine &line, PointType &pos) const
{
    auto &c = line.a;
    auto &d = line.b;
    
    // Line AB represented as a1x + b1y = c1
    double a1 = b.y - a.y;
    double b1 = a.x - b.x;
    double c1 = a1*(a.x) + b1*(a.y);
  
    // Line CD represented as a2x + b2y = c2
    double a2 = d.y - c.y;
    double b2 = c.x - d.x;
    double c2 = a2*(c.x)+ b2*(c.y);
  
    double determinant = a1*b2 - a2*b1;
  
    if (determinant == 0)
    {
        // The lines are parallel. This is simplified
        // by returning a pair of FLT_MAX
        return false;
    }

    double x = (b2*c1 - b1*c2)/determinant;
    double y = (a1*c2 - a2*c1)/determinant;
    pos = PointType(x, y, 0.0f);
    
    return true;
}

void PolyLine::SortIntersectionsList(PolyDetector &pd)
{
    std::sort(intersections.begin(), intersections.end(), [&pd, this] (const uint32_t &p1, const uint32_t &p2) {
        //return bComparePointOrder(pd.intersectionPoints[p1], pd.intersectionPoints[p2]);
        return pd.intersectionPoints[p1].squaredist(a) < pd.intersectionPoints[p2].squaredist(a);
    });
}

PolyLine PolyDetector::newLine(uint32_t i, uint32_t j, PolyLine &origLine)
{
    PolyLine l;
    
    l.aIdx = origLine.intersections[i];
    l.bIdx = origLine.intersections[j];
    l.a = intersectionPoints[l.aIdx];
    l.b = intersectionPoints[l.bIdx];
    l.origLine = origLine.id;
    l.lastDissolveStep = origLine.lastDissolveStep;
    l.attr0 = origLine.attr0;
    l.calcCenter();
    
    return l;
}

bool PolyDetector::CreateLines()
{
    //logoutf("Line creation");
    
    uint32_t n = 0;
    for (auto &l : origLines)
    {
        l.id = n++;
    }
    
    // prior to removing overlapping, one must
    // remove all zero length line, otherwise the results
    // will be unpredictable
    RemoveZeroLengthLines();

    // then we must remove line overlapping
    RemoveOverlappings();

    // finally we detect intersections between lines
    int intersection_count = DetectAllIntersections();
    if (intersection_count == 0)
        return true;
    if (verbose)
    {
        logoutf("Detected %d intersections", intersection_count);
        logoutf("%u lines after intersection detection", uint32_t(origLines.size()));
    }

    // sweep all lines
    lines.clear();
    for (auto &line : origLines)
    {
        if (line.ignore) continue;
        
        // check if current line has intersections
        if (line.intersections.size() >= 2)
        {
            line.SortIntersectionsList(*this);
            if (verbose > 1)
                logoutf("line #%u nIntersections:%u", line.id, uint32_t(line.intersections.size()));
            
            for (uint32_t i = 1; i < line.intersections.size(); ++i)
            {
                bool foundDup = false;
                auto
                    aIdx = line.intersections[i - 1],
                    bIdx = line.intersections[i];
                
                assert(aIdx != bIdx);
                
                auto &p1 = intersectionPoints[aIdx];
                auto &p2 = intersectionPoints[bIdx];
                if (!pointsDiffer(p1, p2))
                {
                    logoutf("P%u P%u are the same for line #%u", aIdx, bIdx, line.id);
                    assert(pointsDiffer(p1, p2));
                }
                
                assert(p1.squaredist(line.a) <= p2.squaredist(line.a));

                for (auto &lDup : lines)
                {
                    if ((lDup.minPid() == std::min(aIdx, bIdx)) &&
                         lDup.maxPid() == std::max(aIdx, bIdx))
                    {
                        //if (verbose)
                            logoutf("WARN: Identical lines detected! procLine:%u [P%u P%u]. Means origLine:#%u and origLine:#%u are collinear!", lDup.id,
                                    lDup.minPid(), lDup.maxPid(), line.id, lDup.origLine);
                        //assert(false);
                        foundDup = true;
                    }
                }
                
                if (foundDup)
                {
                    if (verbose)
                        logoutf("WARN: line #%u has duplicate procLines! [P%u P%u]", line.id, aIdx, bIdx);
                }
                else if (!pointsDiffer(p1, p2))
                {
                    if (verbose)
                        logoutf("WARN: line #%u has zero length! [P%u P%u]", line.id, aIdx, bIdx);
                }
                else
                {
                    lines.push_back(newLine(i - 1, i, line));
                }
                
                if (verbose > 3)
                {
                    logoutf("[L#%u][%f %f, %f %f] -> %u [P%u P%u] [%f %f] [%f %f]",
                            line.id,
                            line.a.x, line.a.y, line.b.x, line.b.y,
                            i,
                            aIdx < bIdx ? aIdx : bIdx, aIdx < bIdx ? bIdx : aIdx,
                            p1.x, p1.y,
                            p2.x, p2.y);
                }
            }
        }
    }
    
#if 1
    for (uint32_t i = 0; i < lines.size(); i++)
    {
        auto &l1 = lines[i];
        
        //if (!pointsDiffer(l1.a, l1.b)) // find a zero length line
        //    logoutf("line %s points are the same!", l1.toString(*this).c_str());
        
        for (uint32_t j = i + 1; j < lines.size(); j++)
        {
            auto &l2 = lines[j];
            
            auto aIdx1 = l1.minPid();
            auto bIdx1 = l1.maxPid();
            
            auto aIdx2 = l2.minPid();
            auto bIdx2 = l2.maxPid();
            
            if (!l1.HasCommonIdxPoints(l2) && overlap(l1, l2))
            {
                logoutf("%u:%s and %u:%s overlap! olid:#%u olid:#%u", i, l1.toString(*this).c_str(), j, l2.toString(*this).c_str(), l1.origLine, l2.origLine);
            }
            
            if (aIdx1 == aIdx2 && bIdx1 == bIdx2)
            {
                logoutf("%u:%s and %u:%s have same points! [P%u P%u] olid1:#%u olid2:#%u", i, l1.toString(*this).c_str(), j, l2.toString(*this).c_str(), aIdx1, bIdx1, l1.origLine, l2.origLine);
                //assert(0);
                
                //dumpLines("commonPts");
            }
        }
    }
#endif
    
    if (verbose)
        logoutf("nOrigLines:%u nLines:%u", uint32_t(origLines.size()), uint32_t(lines.size()));

    return true;
}

void PolyDetector::dumpLines(const char *msg, bool useIgnore)
{
    for (auto &l : lines)
    {
        if (useIgnore && l.ignore) continue;
        logoutf("[%s] %s", msg, l.toString(*this).c_str());
    }
}

/***
* @desc removes all lines with zero length
*/
void PolyDetector::RemoveZeroLengthLines(void)
{
    for (auto it = origLines.begin(); it != origLines.end(); )
    {
        if (!pointsDiffer(it->a, it->b)) // find a zero length line
            it = origLines.erase(it);
        else
            ++it;
    }
}

bool PolyDetector::addPointToLine(uint32_t pid, uint32_t lid)
{
    auto &v = pointToLines[pid];
    if (std::find(v.begin(), v.end(), lid) == v.end())
    {
        v.push_back(lid);
        return true;
    }
    return false;
}

/***
* @descr removes line overlappings
* @note must be called before applying Bentley-Ottmann algorithm
*/
void PolyDetector::RemoveOverlappings()
{
    uint32_t i, j, count = uint32_t(origLines.size());
    
    //logoutf("RemoveOverlappings");
    
    for (auto &l : origLines)
    {
        l.calcCenter();
    }
    
    PolyLine line;

    // lets find overlapping lines
    uint32_t countBefore = count;
    for (i = 0; i < count; i++)
    {
        auto &line_i = origLines[i];

        for(j = i + 1; j < count; j++)
        {
            auto &line_j = origLines[j];
            if (overlap(line_i, line_j))
            {
                int ret = simplifiedLine(line_i, line_j, line);
                
                if (verbose > 3)
                    logoutf("origLine #%u {%f %f %f %f} overlaps with origLine #%u {%f %f %f %f} ret:%d",
                            line_i.id, line_i.a.x, line_i.a.y, line_i.b.x, line_i.b.y,
                            line_j.id, line_j.a.x, line_j.a.y, line_j.b.x, line_j.b.y,
                            ret);

                if (ret == 1)
                {
                    // must remove line_j
                    if (verbose > 1)
                        logoutf("rm origLine %u", line_j.id);
                    origLines.erase(origLines.begin() + j);
                    j--;
                    count--;
                    
                    //i = 0; break;
                }
                else
                {
                    if (ret != 2)
                    {
                        //line.id = uint32_t(origLines.size());
                        line.id = line_i.id;
                        
                        if (verbose)
                            logoutf("new origLine %u by merging %u and %u", line.id, line_i.id, line_j.id);
                        
                        // must remove both line_i and line_j and add a new one
                        origLines.erase(origLines.begin() + j);
                        
                        line.calcCenter();
                        origLines.push_back(std::move(line));
                    }

                    // must remove line_i
                    origLines.erase(origLines.begin() + i);

                    // update counters
                    i--;
                    count--;

                    // skip inner loop an go to next step of outer loop
                    break;
                }
            }
        }
    }
    
    if (countBefore != count)
        logoutf("orig: countBefore:%u count:%u", countBefore, count);
}

uint32_t PolyDetector::DetectAllIntersections()
{
    uint32_t ret = 0;
    size_t counter = origLines.size();
    PointType intersection;
    
    intersectionPoints.clear();
    collinearLineMap.clear();
    
    for (auto &l : origLines)
    {
        l.calcCenter();
    }
    
    // sort to always have the same results
    std::sort(origLines.begin(), origLines.end(), PolyLine::bCompareLineOrder);
    
    uint32_t n = 0;
    for (auto &l : origLines)
    {
        l.id = n++;
        l.intersections.clear();
        l.intersectedLines.clear();
    }
    
    // intersected lines: remove lines with only one intersection
    for (uint32_t i = 0; i < counter; i++)
    {
        auto &l1 = origLines[i];
        for (uint32_t j = i + 1; j < counter; j++)
        {
            auto &l2 = origLines[j];
            if (doIntersect(l1.a, l1.b, l2.a, l2.b))
            {
                //if (verbose > 1)
                //    logoutf("line #%u intersects #%u", l1.id, l2.id);
                
                /*
                for (auto &l : {&l1, &l2})
                {
                    auto &otherLid = l->id == l1.id ? l2.id : l1.id;

                    if (l->intersectedLines.find(otherLid) != l->intersectedLines.end())
                    {
                        logoutf("WARN: origLine #%u already contains intersected origLine #%u", l->id, otherLid);
                    }
                    else
                    {
                        l->intersectedLines.insert(otherLid);
                    }
                }
                */
                
                //if (l1.intersectedLines.find(l2.id) != l1.intersectedLines.end())
                l1.intersectedLines.insert(l2.id);
                l2.intersectedLines.insert(l1.id);
            }
        }
    }
    
    for (uint32_t i = 0; i < counter; i++)
    {
        auto &l1 = origLines[i];
        if (l1.intersectedLines.size() == 1)
        {
            if (verbose)
                logoutf("line #%u has only %u intersectedLines", l1.id, uint32_t(l1.intersectedLines.size()));
            
            for (uint32_t j = 0; j < counter; j++)
            {
                if (i != j)
                {
                    auto &l2 = origLines[j];
                    
                    for (auto it = l2.intersectedLines.begin(); it != l2.intersectedLines.end(); )
                    {
                        if (*it == l1.id)
                        {
                            if (verbose)
                                logoutf("removed intersectedLine line #%u from line #%u", l1.id, l2.id);
                            it = l2.intersectedLines.erase(it);
                            i = 0; // recheck all lines
                        }
                        else
                            ++it;
                    }
                }
            }
            
            l1.intersectedLines.clear();
        }
    }
    
    if (verbose > 2)
    {
        for (auto &l : origLines)
        {
            //if (!l.intersectedLines.empty())
            {
                std::string str;
                for (auto &n : l.intersectedLines)
                    str += std::to_string(n) + " ";
                if (!str.empty())
                    str.pop_back();
                logoutf("line #%u has %u intersectedLines: [%s]", l.id, uint32_t(l.intersectedLines.size()), str.c_str());
            }
        }
    }
    
    // intersection points
    std::set<std::pair<uint32_t, uint32_t>> took;
    for (auto &l1 : origLines)
    {
        for (auto &lid2 : l1.intersectedLines)
        {
            auto &l2 = origLines[lid2];
            
            if (took.find(std::make_pair(l1.id, l2.id)) == took.end())
            {
                intersection = vec(0);
                if (l1.IntersectionPoint(l2, intersection)) // checks if not parallel
                {
                    uint32_t intersectionIdx = uint32_t(intersectionPoints.size());
                    
                    // Check if the same intersection point exists
                    bool dupPoint = false;
                    for (uint32_t pi = 0; pi < intersectionPoints.size(); ++pi)
                    {
                        auto &p = intersectionPoints[pi];
                        if (!pointsDiffer(p, intersection))
                        {
                            if (verbose > 3)
                            {
                                logoutf("WARN: origLine #%u intersects origLine #%u exactly in the same point P%u as:", l1.id, l2.id, pi);
                                for (auto &lDup : origLines)
                                {
                                    if (lDup.id != l1.id && lDup.id != l2.id)
                                    {
                                        for (auto &inters : lDup.intersections)
                                        {
                                            if (inters == pi)
                                            {
                                                logoutf("WARN:  origLine #%u", lDup.id);
                                            }
                                        }
                                    }
                                }
                            }
                            dupPoint = true;
                            intersectionIdx = pi;
                            intersection = p;
                            break;
                        }
                    }

                    if (!dupPoint)
                    {
                        intersectionPoints.push_back(intersection);
                    }
                    
                    if (verbose > 2)
                        logoutf("[P%u]: origLine #%u intersects origLine #%u. dup:%d", intersectionIdx, l1.id, l2.id, dupPoint);
                    
                    for (auto &l : {&l1, &l2})
                    {
                        if (std::find(l->intersections.begin(), l->intersections.end(), intersectionIdx) != l->intersections.end())
                        {
                            if (verbose > 1)
                                logoutf("WARN: duplicate intersectionIdx P%u in line #%u", intersectionIdx, l->id);
                        }
                        else
                        {
                            //auto d = obb::LineSegment(l->a, l->b).Distance(intersection);
                            //assert(d < .1f);
                            
                            //auto d1 = obb::LineSegment(l->a, l->b).Distance(intersectionPoints[intersectionIdx]);
                            //assert(d1 < .1f);
                            
                            l->intersections.push_back(intersectionIdx);
                        }
                    }
                    
                    //addPointToLine(intersectionIdx, l1.id);
                    //addPointToLine(intersectionIdx, l2.id);
                    
                    took.insert(std::make_pair(l1.id, l2.id));
                    took.insert(std::make_pair(l2.id, l1.id));
                    
                    ret++;
                }
            }
        }
    }
    if (verbose)
        logoutf("intersectionPoints:%u", uint32_t(intersectionPoints.size()));
    
    //return 0;
    
    if (verbose > 2)
    {
        for (auto &l : origLines)
        {
            //if (!l.intersections.empty())
            {
                std::string str;
                for (auto &n : l.intersections)
                    str += std::to_string(n) + " ";
                if (!str.empty())
                    str.pop_back();
                logoutf("line #%u has %u pts: [%s]", l.id, uint32_t(l.intersections.size()), str.c_str());
            }
        }
    }
    
#if 1
    
    uint32_t nCol = 0;
    std::vector<uint32_t> pids;
    
    uint32_t times = 0;
    bool ok = false;
    do
    {
        took.clear();
        //nCol = 0;
        
        ok = true;
        
        for (auto &l1 : origLines)
        {
            if (l1.ignore) continue;
            if (l1.intersections.size() < 2) continue;
            
            //if (times == 1)
            //    l1.SortIntersectionsList(*this);
            
            float a, b, c;
            vec::line(vec(l1.a.x, l1.a.y), vec(l1.b.x, l1.b.y), a, b, c);
            if (a + b == 0.0f)
            {
                logoutf("l1:#%u l1.a: [%f %f] l1.b: [%f %f]", l1.id, l1.a.x, l1.a.y, l1.b.x, l1.b.y);
                assert(0);
                continue;
            }
            
            //obb::Line line = obb::LineSegment(l1.a, l1.b).ToLine();
        
            for (auto &l2 : origLines)
            {
                if (l2.ignore) continue;
                if (l2.intersections.size() < 2) continue;
                
                //if (took.find(std::make_pair(l1.id, l2.id)) == took.end())
                {
                    //if (k == 1)
                    //    l2.SortIntersectionsList(*this);
                    
                    took.insert(std::make_pair(l1.id, l2.id));
                    took.insert(std::make_pair(l2.id, l1.id));
                    
                    if (l1.id == l2.id) continue;
                    if (l2.intersections.empty()) continue;
                    assert(&l1 != &l2);
                    
                    uint32_t nFound = 0;
                    float maxLineDist = 0;
                    pids.clear();
                    for (auto &pid1 : l1.intersections)
                    {
                        for (auto &pid2 : l2.intersections)
                        {
                            if (pid1 == pid2)
                            {
                                nFound++;
                                pids.push_back(pid1);
                            }
                            
                            //float d = line.Distance(intersectionPoints[pid2]);
                            float d = vec::lineDist(a, b, c, vec(intersectionPoints[pid2].x, intersectionPoints[pid2].y));
                            if (d > maxLineDist)
                                maxLineDist = d;
                        }
                    }
                    
                    if (nFound >= 2)
                    {
                        ok = false;
                        
#if 1
                        if (maxLineDist <= minPointDiff)
                        {
                            logoutf("[%u]: line #%u is collinear with line #%u! nFound:%u lineDist:%f merging ...", times, l1.id, l2.id, nFound, maxLineDist);
                            
                            // TODO: merge point of l2 into l1
                            uint32_t nMerged = 0;
                            for (auto &pid : l2.intersections)
                            {
                                //logoutf("add P%u to l#%u", pid, l1.id);
                                if (std::find(l1.intersections.begin(), l1.intersections.end(), pid) == l1.intersections.end())
                                {
                                    //logoutf("add P%u to l#%u", pid, l1.id);
                                    l1.intersections.push_back(pid);
                                    nMerged++;
                                }
                            }
                            if (nMerged > 0)
                            {
                                logoutf("[%u]: merged %u points from line #%u into line #%u! dist:%f", times, nMerged, l2.id, l1.id, maxLineDist);
                            }
                            
                            //logoutf("line #%u with %u intersections ignored", l2.id, uint32_t(l2.intersections.size()));
                            l2.intersections.clear();
                            l2.a = l2.b = l2.center = vec(0);
                            l2.ignore = true;
                            
                            nCol++;
                        }
                        else
                        {
                            // merge points?
                            
                            // remove them
                            uint32_t nRm1 = 0, nRm2 = 0;
                            for (auto &pid : pids)
                            {
                                logoutf("[%u]]: P%u", times, pid);
                                for (auto it = l1.intersections.begin(); it != l1.intersections.end();)
                                {
                                    if (*it == pid)
                                    {
                                        it = l1.intersections.erase(it);
                                        nRm1++;
                                    }
                                    else
                                        ++it;
                                }
                                for (auto it = l2.intersections.begin(); it != l2.intersections.end();)
                                {
                                    if (*it == pid)
                                    {
                                        it = l2.intersections.erase(it);
                                        nRm2++;
                                    }
                                    else
                                        ++it;
                                }
                            }
                            auto keepPid = pids[0];
                            l1.intersections.push_back(keepPid);
                            l2.intersections.push_back(keepPid);
                            
                            logoutf("[%u]]: line #%u has common points(%u) with line #%u keeping just P%u nRm1:%u nRm2:%u rem1:%u rem2:%u lineDist:%f",
                                    times, l1.id, uint32_t(pids.size()), l2.id, keepPid, nRm1, nRm2,
                                    uint32_t(l1.intersections.size()),
                                    uint32_t(l2.intersections.size()),
                                    maxLineDist);
                        }
#endif
                    
                    }
                }
            }
        }
        if (ok)
            logoutf("[%u] ok!", times);
        times++;
    }
    while (!ok);
    
    logoutf("Number of col lines using intersection points: %u", uint32_t(nCol));
#endif

    for (auto &l : origLines)
    {
        if (l.intersections.size() == 1)
        {
            logoutf("WARN: line #%u has only %u intersections!", l.id, uint32_t(l.intersections.size()));
            //assert(false);
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
        
        std::swap(aIdx, bIdx);
    }
}

/***
* @descr sort the lines
*/
void PolyDetector::SortLines(void)
{
    for (auto &line : lines)
    {
        if (!line.ignore)
            line.CalculateFirstAndLastPoint();
    }
    std::sort(lines.begin(), lines.end(), PolyLine::bCompareLineOrder);
}

std::string PolyLine::toString(PolyDetector &pd) const
{
    uint32_t nNeigh = 0;
    std::string neighStr = neighToString(pd, &nNeigh);
    std::string str = "L{[" + std::to_string(id) + "] nNeigh:";
    str += std::to_string(nNeigh);
    
    str += " [";
    str += neighToString(pd);
    str += "] [P";
    str += std::to_string(aIdx < bIdx ? aIdx : bIdx);
    str += " P";
    str += std::to_string(aIdx < bIdx ? bIdx : aIdx);
    str += ']';
    
    str += " took:" + std::to_string(took);
    
    if (ignore)
        str += " ING!";
    
    str += '}';
    
    return str;
}

std::string PolyLine::neighToString(PolyDetector &pd, uint32_t *retNNeigh) const
{
    auto search = pd._neighbors.find(id);
    if (search == pd._neighbors.end())
        return "[]";
    std::string str;
    uint32_t nNeigh = 0;
    for (auto &n : search->second)
    {
        auto l = pd.findLine(n);
        if (l)
        {
            str += std::to_string(n) + " ";
            nNeigh++;
        }
    }
    if (!str.empty())
        str.pop_back();
    if (retNNeigh)
        *retNNeigh = nNeigh;
    return str;
}

uint32_t PolyLine::numNeigh(PolyDetector &pd) const
{
    auto search = pd._neighbors.find(id);
    if (search == pd._neighbors.end())
        return 0;
    uint32_t ret = 0;
    for (auto &n : search->second)
    {
        auto l = pd.findLine(n);
        if (l)
            ret++;
    }
    return ret;
}
uint32_t PolyLine::numIntersections(PolyDetector &pd) const
{
    return uint32_t(intersections.size());
}

uint32_t PolyLine::canBeRemoved(PolyDetector &pd, RmLinesType type) const
{
    assert(!ignore);
    
    //if (type == 1 && verbose)
    //    logoutf("canBeRemoved(%s) type:%d", toString(pd).c_str(), type);
    
    if (type == RmLinesType::NoPointNeigh)
    {
        // Remove lines having points with no neighbors
        for (auto &abId : {aIdx, bIdx})
        {
            auto &neigh = pd.pointToLines[abId];
            if (neigh.size() <= 1) // including self
            {
                //if (pd.verbose)
                //    logoutf("canBeRemoved(%s) point with no neighbors! CAN be removed", toString(pd).c_str());
                return true;
            }
            
            // check also ignore flag
            uint32_t nNeigh = 0;
            for (auto &n : neigh)
            {
                auto l = pd.findLine(n);
                if (l && l->id != id)
                    nNeigh++;
            }
            if (nNeigh == 0)
            {
                if (pd.verbose)
                    logoutf("canBeRemoved(%s) point P%u with no neighbors (%u) (check ignore flag)! CAN be removed", toString(pd).c_str(), abId, nNeigh);
                return true;
            }
        }
    }
    if (type == RmLinesType::TakenTwice)
    {
        if (took >= 2)
        {
            for (auto &lid : pd._neighbors[id])
            {
                auto l = pd.findLine(lid);
                if (l && l->id != id)
                {
                    //if (pd.verbose)
                    //    logoutf("nlid %s", l->toString(pd).c_str());
                    
                    //if (!l->took)
                    if (l->took < 0)
                    {
                        return false;
                    }
                }
            }
            //logoutf("canBeRemoved(%s) line tooked 2 times! CAN be removed", toString(pd).c_str());
            return true;
        }
    }
    
    return false;
}

bool PolyDetector::DetectPolygons()
{
    polys.clear();
    
    logoutf("verbose:%u", verbose);
    
    if (verbose)
        logoutf("origLines: %u", uint32_t(origLines.size()));

    if (!CreateLines())
    {
        logoutf1("Error: Could not successfully create lines!");
        return false;
    }
    if (verbose)
        logoutf("lines %u", uint32_t(lines.size()));
    
    while (true)
    {
        bool stop = false;
        
        if (verbose)
            logoutf("----------- BEGIN %u {", dissolveCount);
        
        SortLines();
        
        if (!FindPolys())
        {
            logoutf1("Error constructing the polygon set.");
            return false;
        }

        if (verbose)
            logoutf("Polygon set contains %d polygons.", GetPolyCount());
        
        SimplifyPolys(0.0);
        //logoutf("Polygon set contains %d polygons after simplification!.", GetPolyCount());
        
        auto beforeCount = dissolveCount;
        
        if (verbose)
            logoutf("----------- DISSOLVE %u", dissolveCount);
        
        if (!dissolve())
            stop = true;
        
        if (verbose)
            logoutf("} ----------- END %u\n", beforeCount);
        
        if (stop)
            break;
    }
    
#if 1
    //std::sort(polys.begin(), polys.end(), [](const PolyPol &a, const PolyPol &b) {
    //    return a.dissolveStep < b.dissolveStep;
    //});
    
    std::string str;
    
    for (auto &p : polys)
    {
        //if (verbose)
        //    p.cycle.print(("[FINI_ALL STEP:" + std::to_string(p.dissolveStep) + "] ").c_str());
        
        /*
        str = "[";
        for (auto &lid : p.cycle.idx)
        {
            auto l = findLine(lid);
            if (!l) continue;
            
            str += "{";
            for (auto &pid : {l->aIdx, l->bIdx})
                str += std::to_string(pid) + " ";
            str.pop_back();
            str += "} ";
        }
        if (!str.empty())
            str.pop_back();
        str += "]";
        
        logoutf("[%u]: pts:%s", p.dissolveStep, str.c_str());
        */
    }
#endif

    return true;
}

/*
bool PolyPol::IsClosed() const
{
    return !p.empty() && !pointsDiffer(p[0], p.back());
}
*/

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
/*
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
*/

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
    
    /*
    if (!IsClosed())
    {
        p.clear();
        return;
    }
    */
    
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

/***
* @desc simplifies the polygons in this set
* @note removes inclusions and disposes small polygons
*/
void PolyDetector::SimplifyPolys(double smaller_polygon_length)
{
    //logoutf("Polygon set simplification");

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

    if (nRemoved)
        logoutf("Merged %d contained polygons", nRemoved);
#endif
}

void PolyDetector::AddLine(const PolyLine &line)
{
    //logoutf("add line #%u (%f %f %f %f)", uint32_t(origLines.size()), line.a.x, line.a.y, line.b.x, line.b.y);
    origLines.push_back(line);
}

void PolyPol::addLine(const PolyLine &l)
{
    p.push_back(l.a);
    p.push_back(l.b);
}

PolyLine *PolyDetector::findLine(uint32_t pidA, uint32_t pidB, bool useIgnore)
{
    auto m = std::min(pidA, pidB);
    auto M = std::max(pidA, pidB);
    
#if 1
    for (auto &l : lines)
    {
#else // faster using neighbors, but they need to be populated
    for (auto &lid : _neighbors[pidA])
    {
        auto lp = findLine(lid, useIgnore);
        if (!lp) continue;
        auto &l = *lp;
#endif
        if (useIgnore && l.ignore) continue;
        if (l.minPid() == m && l.maxPid() == M)
            return &l;
    }
    return nullptr;
}

PolyLine *PolyDetector::findLine(uint32_t id, bool useIgnore)
{
    auto search = lineIdToIdx.find(id);
    if (search != lineIdToIdx.end())
    {
        //assert(search->second >= 0 && search->second < lines.size());
        auto &l = lines[search->second];
        if (useIgnore && l.ignore)
        {
            //logoutf("line id %u is ignored!", id);
            return nullptr;
        }
        
        return &l;
    }
    if (search == lineIdToIdx.end())
    {
        logoutf("WARN: line %u can't be found in lineIdToIdx! szFast:%u szLines:%u v:%d", id,
                uint32_t(lineIdToIdx.size()),
                uint32_t(lines.size()),
                search != lineIdToIdx.end() ? search->second : -2);
    }
    
    return nullptr;
}

PolyLine *PolyDetector::findOrigLine(uint32_t id)
{
    for (auto &l : origLines)
    {
        if (l.id == id)
        {
            return &l;
        }
    }
    logoutf("Cannot find origLine with id %u. nlines:%u", id, uint32_t(origLines.size()));
    return nullptr;
}

bool similarCycle(const PolyCycles &cycles, const PolyCycle &cycle)
{
#if 1
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
#else
    if (cycle.idx.size() < 3) return true;
    
    for (auto &c : cycles)
    {
        assert(c.idx.size() >= 3);
        auto end1 = c.idx.cend();
        auto end2 = cycle.idx.cend();
        --end1;
        --end2;
        if (*c.idx.cbegin() == *cycle.idx.cbegin() && *end1 == *end2)
            return true;
        
        if (c.idx.size() > cycle.idx.size() && c.contains(*end2) && c.contains(*cycle.idx.cbegin()))
            return true;
        if (c.idx.size() < cycle.idx.size() && cycle.contains(*end1) && c.contains(*c.idx.cbegin()))
            return true;
    }
    return false;
#endif
}

bool PolyCycle::Equals(const PolyCycle &p) const
{
    if (idx.size() != p.idx.size())
        return false;
    for (auto it1 = idx.cbegin(), it2 = p.idx.cbegin(); it1 != idx.cend(); ++it1, ++it2)
    {
        if (*it1 != *it2)
            return false;
    }
    return true;
}

#if 0
uint32_t PolyCycle::numCuts(PolyDetector &pd, const PolyLine &l) const
{
    uint32_t ret = 0;
    for (auto &id1 : idx)
    {
        auto l1 = pd.findLine(id1);
        if (l1)
        {
            if (l.HasCommonIdxPoints(*l1))
                ret++;
        }
    }
    return ret;
}
#endif

bool PolyCycle::AddLineId(PolyDetector &pd, uint32_t id)
{
    auto l = pd.findLine(id);
    if (!l)
    {
        logoutf("findLine(%u) failed!", id);
        return false;
    }
    
    if (l->took >= 2)
    {
        if (pd.verbose > 3)
            logoutf("line %u can't be added to cycle! tooked %u times!", id, l->took);
        return false;
    }
    
    l->test0 = 0; // a cnt
    l->test1 = 0; // b cnt
    
    bool shareA = false, shareB = false;
    bool collinear = false;
    
    for (auto &id1 : idx)
    {
        auto l1 = pd.findLine(id1);
        if (!l1) continue;
        
        shareA = l->aIdx == l1->aIdx || l->aIdx == l1->bIdx;
        shareB = l->bIdx == l1->aIdx || l->bIdx == l1->bIdx;
        collinear = false;

        // share a point
        if (shareA)
            l->test0++;
        
        // share b point
        if (shareB)
            l->test1++;
        
        if (l->test0 >= 2 || l->test1 >= 2)
        {
            if (pd.verbose > 2)
                logoutf("line %u can't be added to cycle! ta:%d tb:%d ids: [%s]", id, l->test0, l->test1, idxToString().c_str());
            return false;
        }
        
        if (shareA || shareB)
        {
            collinear = pd.CollinearIdx(*l, *l1);
            if (pd.verbose > 2)
                logoutf("l(%u), l1(%u) sa:%d sb:%d col:%d", l->id, l1->id, shareA, shareB, collinear);
            if (collinear)
            {
                if (pd.verbose > 2)
                    logoutf("line %u can't be added to cycle! colliniar with %u sa:%d sb:%d", id, id1, shareA, shareB);
                return false;
            }
        }
        
        for (auto &id2 : idx)
        {
            if (id2 != id && id2 != id1
                //&& !pd.CollinearIdx(id2, id) && !pd.CollinearIdx(id2, id1)
            )
            {
                auto l2 = pd.findLine(id2);
                
                if (l2 && doIntersect(l->center, l1->center, l2->a, l2->b))
                //if (l2 && !Collinear(ls, *l2) && ls.PolyIntersects(*l2))
                {
                    if (pd.verbose > 2)
                        logoutf("line %u can't be added to cycle (coll test)! id2:%u intersects with middle(%u, %u)", id, id2, l->id, l1->id);
                    return false;
                }
            }
        }
    }
    
#if 1 // brute force
    for (auto &pid : {l->aIdx, l->bIdx})
    {
        auto &pl = pd.pointToLines[pid];
        if (pl.size() >= 2)
        {
            for (auto it1 = pl.begin(); it1 != pl.end(); ++it1)
            {
                auto &nlid1 = *it1;

                if (nlid1 != id && contains(nlid1))
                {
                    for (auto it2 = pl.begin(); it2 != pl.end(); ++it2)
                    {
                        if (*it1 == *it2) continue;
                        
                        auto nlid2 = *it2;
                        
                        if (nlid2 != id && !contains(nlid2))
                        {
                            //auto l1 = pd.findLine(nlid1, false);
                            auto l1 = pd.findLine(nlid1, false);
                            auto l2 = pd.findLine(nlid2);
                            if (l1 && l2)
                            {
                                //logoutf("checking lines %u %u %u in P%u", l->id, l1->id, l2->id, pid);
                                if (l2->betweenNeighbors(pd, *l1, *l))
                                    return false;
                            }
                        }
                    }
                }
            }
        }
    }
#endif
    
    
    idx.insert(id);
    //l->processed = id;
    lastIdx = id;
    l->processed = startIdx + 1;
    return true;
}

bool PolyCycle::accepted(PolyDetector &pd)
{
    for (auto &lid : idx)
    {
        auto l = pd.findLine(lid);
        if (l)
        {
            l->incTook(pd);
            l->processed = startIdx + 1;
        }
    }
    return true;
}

/*
static bool equalCycles(const CycleSet &a, const CycleSet &b)
{
    if (a.size() != b.size())
        return false;

#if 0
    return a == b;
#elif 0
    

#if 1
    return std::equal(a.begin(), a.end(), b.begin());
#else
    for (auto aIt = a.begin(), bIt = b.begin(); aIt != a.end(); ++aIt, ++bIt)
        if (*aIt != *bIt)
            return false;
    return true;
#endif
    
#else
    
    assert(a.size() >= 3);
    assert(b.size() >= 3);
    
    auto end1 = a.cend();
    auto end2 = b.cend();
    --end1;
    --end2;
    if (*a.cbegin() == *b.cbegin() && *end1 == *end2)
        return true;

    return false;

#endif // ==
}
*/

/*
bool PolyDetector::cycleProcessed(const PolyCycle &cycle) const
{
    if (cycle.idx.size() < 3)
        return false;
    
    for (auto &c : processed)
        if (c.idx.size() == cycle.idx.size() && equalCycles(c.idx, cycle.idx))
            return true;
    return false;
}
 */

bool PolyDetector::BuildCycle(uint32_t id, PolyCycle cycle) // as value!
{
    //if (cycle.startIdx != id && cycle.idx.size() > 2)
    //    return true;
    
    //if (cycle.idx.size() >= 7)
    //    return true;
    
    auto l = findLine(id);
    if (!l)
        return true;
    
    if (_neighbors[id].size() < 2)
    {
        return true;
    }
    
    if (verbose > 2)
        cycle.print((std::string("[PROC:") + std::to_string(id) + std::string("]")).c_str());
    
    if (cycle.canBeClosed(*this, id))
    {
        cycle.isClosed = true;
        
        if (verbose > 2)
            cycle.print("[CLOSED] ");
        
        if (!similarCycle(_cycles, cycle))
        {
            if (cycle.accepted(*this))
            {
                _cycles.push_back(cycle);
                if (verbose > 2)
                    cycle.print("[ACCEPTED] ");
                //logoutf("nCycles:%u", uint32_t(_cycles.size()));
                //return false;
            }
            else
            {
                if (verbose > 2)
                    cycle.print("[NOT ACCEPTED] ");
            }
        }
        else
        {
            if (verbose > 2)
                cycle.print("[CYCLE_EXISTS] ");
        }
        
        return true;
    }
    
    if (!cycle.AddLineId(*this, l->id))
    {
        return true;
    }
    
#if 0 // disable entire process checking
#if 0
    if (l->processed + 1 == cycle.startIdx)
    {
        if (verbose > 2)
        {
            logoutf("line %s can't be added to cycle %s! processed: %u!", l->toString(*this).c_str(), cycle.toString().c_str(), l->processed - 1);
        }
        //printf("procHit for %u!\n", l->id);
        return true;
    }
#else
    if (cycleProcessed(cycle))
    {
        if (verbose > 2)
            cycle.print("[PROCESSED!] ");
        return true;
    }
    processed.push_back(cycle.idx);
#endif
#endif

    for (auto &nid : _neighbors[id])
    {
        if (_neighbors[nid].size() < 2) continue;
        if (cycle.canBeClosed(*this, nid) || !cycle.contains(nid))
        {
            if (verbose > 2)
                cycle.print((std::string("[NEIGN:") + std::to_string(nid) + std::string("]")).c_str());
            if (!BuildCycle(nid, cycle))
            {
                return false;
                //break;
            }
        }
    }
    
    return true;
}

bool PolyLine::betweenNeighbors(PolyDetector &pd, const PolyLine &l1, const PolyLine &l2) const
{
    auto cpid = commonPid(l1);
    if (cpid < 0)
    {
        //logoutf("l1: cpid %d between %u and %u", cpid, id, l1.id);
        //logoutf("%s %s", toString(pd).c_str(), l1.toString(pd).c_str());
        return false;
    }
    auto cpid2 = commonPid(l2);
    if (cpid != cpid2)
    {
        //logoutf("l2: cpid %d between %u and %u", cpid2, id, l2.id);
        return false;
    }
    
    //auto &p1 =  pd.intersectionPoints[l1.otherPid(cpid)];
    //auto &p2 =  pd.intersectionPoints[l2.otherPid(cpid)];
    
    vec cp = pd.intersectionPoints[cpid];
    vec p = pd.intersectionPoints[otherPid(cpid)];

    //bool ret = obb::Line(cp, cp.sub(p).normalize()).Intersects(obb::LineSegment(l1.center, l2.center));
    bool ret = doIntersect(cp, p, l1.center, l2.center);
    
    if (ret)
    {
        if (pd.verbose > 2)
        {
            //logoutf("line %s is between %s and %s (commonP: P%d)", toString(pd).c_str(), l1.toString(pd).c_str(), l2.toString(pd).c_str(), cpid);
            logoutf("line %u is between %u and %u (commonP: P%d)", id, l1.id, l2.id, cpid);
        }
    }
    else
    {
        //logoutf("line %u is NOT between %u and %u (commonP: P%d)", id, l1.id, l2.id, cpid);
    }
    
    return ret;
}

/*
float PolyLine::angle(PolyDetector &pd, const PolyLine &l) const
{
    auto cpid = l.commonPid(*this);
    
    auto pid = l.otherPid(cpid);
    
    const vec &cp = pd.intersectionPoints[cpid];
    const vec &p = pd.intersectionPoints[pid];
    
    vec dir = vec(b).sub(a);
    vec dir1 = vec(p).sub(cp);
    
    dir.normalize();
    dir1.normalize();
    
    //return dir1.dot(dir);
    return acosf(dir1.dot(dir));
}
*/

bool PolyLine::compareNeigh(PolyDetector &pd, uint32_t nid1, uint32_t nid2) const
{
    auto nl1 = pd.findLine(nid1);
    auto nl2 = pd.findLine(nid2);
    if (!nl1 || !nl2)
    {
        assert(false);
        return true;
    }
    
#if 1
    auto dl1 = nl1->center.squaredist(center);
    auto dl2 = nl2->center.squaredist(center);
    return dl1 < dl2;
#elif 0
    return PolyLine::bCompareLineOrder(*nl1, *nl2);
#else
    return angle(pd, *nl1) < angle(pd, *nl2);
#endif
}

bool PolyLine::sortNeigh(PolyDetector &pd) const
{
    auto &neigh = pd._neighbors[id];
    std::sort(neigh.begin(), neigh.end(), [this, &pd](uint32_t nid1, uint32_t nid2)
    {
        return compareNeigh(pd, nid1, nid2);
    });
    
    return true;
}

uint32_t &PolyLine::incTook(PolyDetector &pd)
{
    if (took >= 2)
        return took;
    took++;
    if (took > 2)
    {
        logoutf("WARN: line %s taken %u times!", toString(pd).c_str(), took);
    }
    
    return took;
}

bool PolyDetector::FindPolys()
{
    //logoutf("FindPolys");
    
    //processed.clear();
    pointToLines.clear();
    collinearLineMap.clear();
    
    // assign line ids
    uint32_t n = 0;
    lineIdToIdx.clear();
    for (auto &l : lines)
    {
        //logoutf("assign id %u to #%u", n, l.id);
        
        if (dissolveCount == 0) // only on first step
        {
            //logoutf("rename line #%u to %u", l.id, n);
            l.id = n;
        }
        
        if (!l.ignore)
        {
            addPointToLine(l.aIdx, l.id);
            addPointToLine(l.bIdx, l.id);
            
            //logoutf("line %u a:P%u b:P%u len:%f", l.id, l.aIdx, l.bIdx, l.a.dist2(l.b));
        }
        
        /*
        if (l.id != n)
        {
            logoutf("line %u does not match n %u", l.id, n);
            assert(l.id != n);
        }
        */
        lineIdToIdx[l.id] = n; // lineIdToIdx[l.id] = lines[n]
        
        n++;
    }
    
    // build collinearLineMap
    for (auto &lo : origLines)
    {
        for (uint32_t i = 0; i < lines.size(); ++i)
        {
            auto &l1 = lines[i];
            if (!l1.ignore && l1.origLine == lo.id)
            {
                for (uint32_t j = i + 1; j < lines.size(); ++j)
                {
                    auto &l2 = lines[j];
                    if (!l2.ignore && l2.origLine == lo.id)
                    {
                        setCollinear(l1.id, l2.id);
                    }
                }
            }
        }
    }
    
    if (verbose > 1)
    {
        for (auto &kv : collinearLineMap)
        {
            for (auto &c : kv.second)
            {
                logoutf("[%u %u] are collinear", kv.first, c);
            }
        }
    }
    if (verbose)
        logoutf("collinearLineMap size:%u", uint32_t(collinearLineMap.size()));

    // Build neighbors
    std::vector<uint32_t> neigh;
    _neighbors.clear();
    for (uint32_t i = 0; i < lines.size(); ++i)
    {
        auto &l1 = lines[i];
        if (!l1.ignore)
        {
            neigh.clear();
            for (uint32_t j = i + 1; j < lines.size(); ++j)
            {
                auto &l2 = lines[j];
                if (!l2.ignore && l1.HasCommonIdxPoints(l2))
                {
                    neigh.push_back(l2.id);
                }
            }
            if (!neigh.empty())
            {
                //std::sort(neigh.begin(), neigh.end(), [this, &l1](const uint32_t &nid1, const uint32_t &nid2) {
                //    return l1.compareNeigh(*this, nid1, nid2);
                //});
                for (auto &nid : neigh)
                {
                    if (verbose > 1)
                        logoutf("line %u intersects %u", l1.id, nid);
                    
                    _neighbors[l1.id].push_back(nid);
                    _neighbors[nid].push_back(l1.id);
                }
            }
        }
    }
    for (auto &kv : _neighbors)
    {
        auto l = findLine(kv.first);
        if (l)
            l->sortNeigh(*this);
    }
    if (verbose)
        logoutf("neighbors MatSize:%u", uint32_t(_neighbors.size()));
    
    if (1)
    {
        for (uint32_t i = 0; i < lines.size(); ++i)
        {
            auto &l1 = lines[i];
            for (uint32_t j = i + 1; j < lines.size(); ++j)
            {
                auto &l2 = lines[j];
                if (l1.id == l2.id)
                {
                    logoutf("line %s and %s have the same ids!", l1.toString(*this).c_str(), l2.toString(*this).c_str());
                    assert(false);
                }
            }
        }
        
        if (verbose > 2)
        //if (verbose)
        {
            for (auto &l : lines)
            {
                if (l.ignore) continue;
                logoutf("FP: %s", l.toString(*this).c_str());
            }
        }
    }
    
    if (verbose > 1)
    {
        for (auto &kv : _neighbors)
        {
            std::string str;
            for (auto &n : _neighbors[kv.first])
            {
                str += std::to_string(n) + " ";
            }
            logoutf("[RT:%u] [%u] nNeigh:%u [%s]", dissolveCount, kv.first, uint32_t(_neighbors[kv.first].size()), str.c_str());
        }
    }
    
    if (dissolveCount == 0)
    {
        for (auto &l : lines)
            l.processed = 0;
    }
    
    for (auto &kv : _neighbors) // point by point
    {
        if (_neighbors[kv.first].size() < 2) continue;
        
        if (verbose > 2)
            logoutf("----------- BEGIN L:%u {", kv.first);
        
        PolyCycle cycle;
        cycle.startIdx = cycle.lastIdx = kv.first;
        cycle.isClosed = false;
        BuildCycle(kv.first, cycle);
        
        if (verbose > 2)
            logoutf("} ----------- END L:%u", kv.first);
    }
    
    std::sort(_cycles.begin(), _cycles.end(), [](const PolyCycle &a, const PolyCycle &b) {
        return a.idxToString().compare(b.idxToString()) < 0;
    });
    
#if 1
    for (auto &cycle : _cycles)
    {
        cycle.fine = true;
        if (!cycle.convex(*this))
            cycle.fine = false;
    }
#endif
    
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
        
        if (verbose)
            cycle.print("[FINI] ");
        
        if (!cycle.fine) continue;
        
        // other the lines for proper triangulation
        std::vector<uint32_t> lidx;
        lidx.reserve(cycle.idx.size());
        for (auto &lid : cycle.idx)
            lidx.push_back(lid);
        std::sort(lidx.begin(), lidx.end(), [this](const uint32_t &lid1, const uint32_t &lid2) {
            auto l1 = findLine(lid1, false);
            auto l2 = findLine(lid2, false);
            if (!l1 || !l2)
            {
                assert(false);
                return true;
            }
            return PolyLine::bCompareLineOrder(*l1, *l2);
        });
        
        // check if not already exists
        bool exists = false;
        for (auto &p : polys)
        {
            if (p.cycle.Equals(cycle))
            {
                exists = true;
                break;
            }
        }
        if (exists)
        {
            logoutf("WARN: cycle %s already exists!", cycle.toString().c_str());
            continue;
        }
        
        PolyPol poly;
        PointType *last = nullptr;
        for (auto &id : lidx)
        {
            auto lPtr = findLine(id);
            if (lPtr)
            {
                lPtr->CalculateFirstAndLastPoint();
                
                auto &a = intersectionPoints[lPtr->aIdx];
                auto &b = intersectionPoints[lPtr->bIdx];
                
                //lPtr->incTook(*this);
                
                poly.addPointChecked(a);
                
                last = &b;
            }
        }
        if (last)
            poly.addPointChecked(*last);
        
        poly.id = cycle.startIdx;
        poly.cycle = cycle;
        
        poly._area = poly.TriangleArea(*this);
        poly.dissolveStep = dissolveCount;
        
        polys.push_back(std::move(poly));
    }
    
    return true;
}

/*
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
*/

/*
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
                    auto pt = l->ClosestPoint(*l1);
                    pt.y -= .01f;
                    p.push_back(pt);
                    break;
                }
            }
        }
    }
}
*/

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

bool PolyPol::addPointChecked(const PointType &v)
{
    for (auto &pt : p)
        if (!pointsDiffer(pt, v))
            return false;
    p.push_back(v);
    return true;
}

/*
obb::Polygon &PolyPol::updatePoly()
{
    poly.p.clear();
    for (auto &pt : p)
        poly.p.push_back(pt);
    return poly;
}
*/

void PolyLine::calcCenter()
{
    CalculateFirstAndLastPoint();
    center = a;
    center.add(b);
    center.mul(0.5f);
}

bool PolyDetector::dissolveCollinear(PolyLine &l1, PolyLine &l2)
{
    if (verbose)
        logoutf("Dissolving collinear neighbors %s and %s", l1.toString(*this).c_str(), l2.toString(*this).c_str());
    
    assert(&l1 != &l2);
    
    uint32_t *p1 = nullptr, *p2 = nullptr;
    if (l1.aIdx == l2.aIdx || l1.aIdx == l2.bIdx)
    {
        p1 = &l1.bIdx; // other point (B)
        p2 = l1.aIdx == l2.aIdx ? &l2.bIdx : &l2.aIdx; // aidx is common point. other point
    }
    else if (l1.bIdx == l2.aIdx || l1.bIdx == l2.bIdx)
    {
        p1 = &l1.aIdx; // other point (A)
        p2 = l1.bIdx == l2.aIdx ? &l2.bIdx : &l2.aIdx; // bidx is common point. other point
    }
    if (!p1 || !p2)
    {
        assert(p1 && p2);
        return false;
    }

#if 1
    auto search = findLine(*p1, *p2, false);
    if (search)
    {
        logoutf("A line %s origLine:#%u already exists with P%u P%u!", search->toString(*this).c_str(), search->origLine, *p1, *p2);
        //assert(false);
        
        l1.setIgnore(*this, "rmCollinear.l1 (dup)");
        l2.setIgnore(*this, "rmCollinear.l2 (dup)");
        
        return false;
    }
#endif

    PolyLine nl;
    nl.id = uint32_t(lines.size());
    nl.aIdx = *p1;
    nl.bIdx = *p2;
    nl.a = intersectionPoints[nl.aIdx];
    nl.b = intersectionPoints[nl.bIdx];
    nl.calcCenter(); // calls CalculateFirstAndLastPoint
    nl.origLine = l1.origLine; // will be used to construct collinearLineMap
    nl.lastDissolveStep = l1.lastDissolveStep;
    nl.attr0 = l1.attr0;
    
    // build point to lines link
    addPointToLine(nl.aIdx, nl.id);
    addPointToLine(nl.bIdx, nl.id);
    
    // build neighbors
    for (auto &id : {l1.id, l2.id})
    {
        for (auto &n : _neighbors[id])
        {
            if (n != l1.id && n != l2.id)
            {
                auto l = findLine(n);
                if (l)
                {
                    _neighbors[nl.id].push_back(n);
                    _neighbors[n].push_back(nl.id);
                }
            }
        }
    }
    
    // set proper collinearities
    for (auto &kv : collinearLineMap)
    {
        if (kv.first == l1.id || kv.first == l2.id)
        {
            for (auto &c : kv.second)
                setCollinear(nl.id, c);
        }
    }
    
    // after neigh built!
    if (verbose)
        logoutf("ignore line %u and line %u", l1.id, l2.id);
    l1.setIgnore(*this, "rmCollinear.l1");
    l2.setIgnore(*this, "rmCollinear.l2");
    
    //return false; // test
    
    lineIdToIdx[nl.id] = uint32_t(lines.size());
    lines.push_back(nl);
    lines.back().sortNeigh(*this);
    if (verbose)
        logoutf("added new line %s by merging %u and %u", lines.back().toString(*this).c_str(), l1.id, l2.id);
    
    // it can be removed
    //dissolveCollinearLine(lines.back());
    
    return true;
}

bool PolyDetector::dissolveCollinearLine(PolyLine &l)
{
    for (auto id : {l.aIdx, l.bIdx})
    {
        uint32_t nValid = 0;
        PolyLine *l1 = nullptr;
        for (auto &n : pointToLines[id])
        {
            if (n != l.id)
            {
                auto lPtr = findLine(n);
                if (lPtr)
                {
                    if (lPtr->took < 2)
                    //if (lPtr->took == 0)
                    //if (lPtr->numNeigh(*this) < 2)
                    {
                        l1 = lPtr;
                        nValid++;
                    }
                }
            }
        }
        //if (nValid == 1)
        //if (nValid > 0)
        //{
        //    logoutf("CollinearIdx[P%u](%s, %s)) = %d. nValid:%u", id, l.toString(*this).c_str(), l1->toString(*this).c_str(), CollinearIdx(l, *l1), nValid);
        //}
        if (nValid == 1 && CollinearIdx(l, *l1))
        {
            if (verbose)
                logoutf("line %s -> P%u has 2 intersections. otherLine:%u", l.toString(*this).c_str(), id, l1->id);
            return dissolveCollinear(l, *l1);
        }
    }
    return false;
}

bool PolyDetector::dissolve()
{
    //return false; // test
    //if (dissolveCount == 1) return false;
    //if (dissolveCount == 1) return false;
    
    if (verbose > 3)
    {
        for (auto &l : lines)
        {
            if (!l.ignore)
            {
                logoutf("DIS: %s", l.toString(*this).c_str());
            }
        }
    }
    
    uint32_t nLinesBefore = uint32_t(lines.size());
    uint32_t nIgnoredBefore = 0;
    for (auto &l : lines)
        if (l.ignore)
            nIgnoredBefore++;
    
    if (_cycles.empty())
    {
        if (verbose)
            logoutf("No need to dissolve for %u cycles", uint32_t(_cycles.size()));
        return false;
    }
    
    if (verbose)
        logoutf("Dissolve for %u cycles", uint32_t(_cycles.size()));
    
    rmLines(RmLinesType::PointConsumed);
    rmLines(RmLinesType::TakenTwice); // points taken twice
    rmLines(RmLinesType::NoPointNeigh);
    rmLines(RmLinesType::Collinear);
    
#if 0
    return false; // test
#endif
    
    uint32_t nIgnored = 0;
    uint32_t nValid = 0;
    for (auto &l : lines)
    {
        if (!l.ignore)
        {
            nValid++;
        }
        else
        {
            nIgnored++;
        }
    }
    if (verbose)
        logoutf("Remaining lines:%u/%u nIgnored:%u", nValid, nLinesBefore, nIgnored);
    
    if (!nValid)
    {
        if (verbose)
            logoutf1("No valid lines. Skipping ...");
        return false;
    }
    
    if (nIgnored == nIgnoredBefore)
    {
        if (verbose)
            logoutf1("No lines ignored. Skipping ...");
        return false;
    }

#if 0
    // reassign new lines
    LineVector newLines;
    newLines.reserve(lines.size());
    for (auto it = lines.begin(); it != lines.end(); ++it)
    {
        if (!it->ignore)
        {
            newLines.push_back(std::move(*it));
        }
    }
    lines = std::move(newLines);
    
    logoutf("Remaining lines:%u/%u (after filter)", uint32_t(lines.size()), nLinesBefore);
#endif
    
    //logoutf("intersectionPoints:%u", uint32_t(intersectionPoints.size())); // should remain the same
    // also clear prev cycles
    _cycles.clear();
        
#if 0 // test
    polys.clear();
#endif
    
    dissolveCount++;
    
    for (auto &l : lines)
    {
        if (!l.ignore)
            l.lastDissolveStep = dissolveCount;
        //l.took = 0;
        if (l.took > 2)
        {
            logoutf("WARN: line %s tooked %u times!", l.toString(*this).c_str(), l.took);
            //assert(false);
            //l.took = 2;
        }
    }
    
    return true;
}

void PolyDetector::setCollinear(uint32_t l1, uint32_t l2)
{
    collinearLineMap[l1].push_back(l2);
    collinearLineMap[l2].push_back(l1);
}

bool PolyDetector::CollinearIdx(uint32_t l1, uint32_t l2)
{
    auto search = collinearLineMap.find(l1);
    if (search == collinearLineMap.end())
        return false;
    return std::find(search->second.begin(), search->second.end(), l2) != search->second.end();
}

bool PolyDetector::CollinearIdx(const PolyLine &l1, const PolyLine &l2)
{
    return CollinearIdx(l1.id, l2.id);
}

bool PolyCycle::canBeClosed(PolyDetector &pd, uint32_t idToAdd) const
{
#if 1
    return idx.size() > 2 && startIdx == idToAdd;
#else
    if (idx.size() < 2) return false;

    auto &neigh = pd._neighbors[startIdx];
    return std::find(neigh.begin(), neigh.end(), idToAdd) != neigh.end();
#endif
}

bool PolyCycle::pointConsumed(PolyDetector &pd, uint32_t pid) const
{
    auto &neigh = pd.pointToLines[pid];
    if (neigh.size() <= 1) // including self
    {
        if (pd.verbose)
            logoutf("pointConsumed(P%u) point with no neighbors! CAN be removed", pid);
        return true;
    }
    
    uint32_t lid1 = 0, lid2 = 0;
    uint32_t nValid = 0, nTaken = 0;
    for (auto &n : neigh)
    {
        auto l = pd.findLine(n);
        if (l)
        {
            nValid++;
            //if (l->took)
            if (l->took)
            {
#if 0
                //if (collinear)
                // if there is at least one collinear neighbor not taken, fail
                for (auto &nlid : _neighbors[l->id])
                {
                    auto ln = findLine(nlid);
                    //if (ln && ln->took < 2 && CollinearIdx(nlid, l->id))
                    if (ln && ln->took < 2)
                    {
                        return false;
                    }
                }
#endif
                if (nTaken == 0)
                    lid1 = l->id;
                else if (nTaken == 1)
                    lid2 = l->id;
                nTaken++;
            }
        }
    }
    
    //logoutf("P%u has nValid:%u taken:%u", pid, nValid, nTaken);
    
    if (nValid == 2 && nValid == nTaken)
    {
        //logoutf("l1:%u l2:%u", lid1, lid2);
        //if (collinear)
            //*collinear = CollinearIdx(lid1, lid2);
        if (pd.CollinearIdx(lid1, lid2))
            return false;
        
        // must be in the cycle!
        if (std::find(idx.begin(), idx.end(), lid1) == idx.end())
        {
            return false;
        }
        if (std::find(idx.begin(), idx.end(), lid2) == idx.end())
        {
            return false;
        }
        
        return true;
    }
    
    return false;
}

void PolyLine::setIgnore(PolyDetector &pd, const char *msg)
{
    if (ignore)
        return;
    
    if (pd.verbose)
        logoutf("[%s] ignore line %s", msg, toString(pd).c_str());
    
    ignore = true;
}

bool PolyDetector::rmEarPoints()
{
    std::vector<PolyLine *> torm;
    
    std::set<uint32_t> pids;
    for (auto &cycle : _cycles)
    {
        for (auto &lid : cycle.idx)
        {
            auto l = findLine(lid);
            if (l)
            {
                for (auto &pid : {l->aIdx, l->bIdx})
                {
                    if (cycle.pointConsumed(*this, pid))
                    {
                        if (verbose)
                            logoutf("Point P%u is an ear! will remove line %u", pid, l->id);
                        torm.push_back(l);
                    }
                }
            }
        }
    }
    //if (verbose)
    //    logoutf("cyclePids:%u/%u", uint32_t(pids.size()), uint32_t(intersectionPoints.size()));
            
    for (auto &l : torm)
    {
        l->setIgnore(*this, "rmEarPoints");
    }
        
    return true;
}

bool PolyDetector::rmLines(RmLinesType type)
{
    if (verbose)
        logoutf("rmLines type:%s", RmLinesTypeStr(type));
    
    if (type == RmLinesType::Collinear)
    {
        // rm collinear lines
        uint32_t nCollinearRemoved = 0;
        for (uint32_t i = 0; i < lines.size(); ++i)
        {
            auto &l = lines[i];
            if (!l.ignore)
            {
                if (dissolveCollinearLine(l))
                {
                    nCollinearRemoved++;
                    i = 0; // again
                }
            }
        }
        if (verbose)
            logoutf("Collinear removed lines:%u", nCollinearRemoved);
        return true;
    }
    
    if (type == RmLinesType::PointConsumed)
    {
        rmEarPoints();
        return true;
    }

    // else pass type to canBeRemoved
    for (uint32_t i = 0; i < lines.size(); ++i)
    {
        auto &l = lines[i];
        if (!l.ignore && l.canBeRemoved(*this, type))
        {
            if (verbose)
                logoutf("remove line %s (type %d)", l.toString(*this).c_str(), type);
            l.setIgnore(*this, (std::string("rmLines.") + RmLinesTypeStr(type)).c_str());
            
            // reprocess same position
            if (i > 0) i--;
            //i = 0;
        }
    }
    
    return true;
}

bool PolyCycle::convex(PolyDetector &pd) const
{
    for (auto it1 = idx.begin(); it1 != idx.end(); ++it1)
    {
        auto l1 = pd.findLine(*it1, false);
        if (!l1) continue;
        
        auto it2 = it1;
        ++it2;
        
        for (; it2 != idx.end(); ++it2)
        {
            auto l2 = pd.findLine(*it2, false);
            if (!l2) continue;
            
            PolyLine ls(l1->center, l2->center);
            for (auto it3 = idx.begin(); it3 != idx.end(); ++it3)
            {
                if (it3 != it1 && it3 != it2
                    //&& !pd.CollinearIdx(*it3, *it1) && !pd.CollinearIdx(*it3, *it2)
                )
                {
                    auto l3 = pd.findLine(*it3);
                    
                    if (l3 && doIntersect(ls.a, ls.b, l3->a, l3->b))
                    //if (l2 && !Collinear(ls, *l2) && ls.PolyIntersects(*l2))
                    {
                        logoutf("convex(%s): l1:%u intersects l2:%u l3:%u", toString().c_str(), l1->id, l2->id, l3->id);
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

void PolyDetector::reset()
{
    _cycles.clear();
    //processed.clear();
    origLines.clear();
    lines.clear();
    polys.clear();
    _neighbors.clear();
    collinearLineMap.clear();
    intersectionPoints.clear();
    pointToLines.clear();
    dissolveCount = 0;
}
