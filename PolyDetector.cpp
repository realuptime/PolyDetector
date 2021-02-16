#include "PolyDetector.h"

#include <assert.h>

#define arToStr(arg) #arg

const char *RmLinesTypeStr(RmLinesType type)
{
    switch (type) {
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

//static bool CollinearVecs(const PointType &p, const PointType &q, const PointType &r)
//{
//    return orientation(p, q, r) == 0;
//}
  
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

static float Area(const PointType &a, const PointType &b, const PointType &c)
{
    return 0.5 * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
}

bool PolyLine::Equals(const PolyLine &line) const
{
    return
        (!pointsDiffer(line.a, a) || !pointsDiffer(line.a, b)) &&
        (!pointsDiffer(line.b, a) || !pointsDiffer(line.b, b));
}

bool PolyLine::HasCommonIdxPoints(const PolyLine &line) const
{
    return
        aIdx == line.aIdx ||
        aIdx == line.bIdx ||
        bIdx == line.aIdx ||
        bIdx == line.bIdx;
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
    return doIntersect(a, b, line.a, line.b);
}

bool PolyLine::IntersectionPoint(const PolyLine &line, PointType &pos) const
{
    return LineLineIntersectionPoint(line, pos);
}

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
    pos =  PointType(x, y, 0.0f);
    
    return true;
}

void PolyLine::SortIntersectionsList(PolyDetector &pd)
{
    std::sort(intersections.begin(), intersections.end(), [&pd] (const uint32_t &p1, const uint32_t &p2) {
        return bComparePointOrder(pd.intersectionPoints[p1], pd.intersectionPoints[p2]);
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
    //logoutf1("Line creation");
    
    // prior to removing overlapping, one must
    // remove all zero length line, otherwise the results
    // will be unpredictable
    RemoveZeroLengthLines();

    // then we must remove line overlapping in order to run
    // the Bentley-Ottmann Algorithm
    RemoveOverlappings();

    // finally we detect intersections between lines
    int intersection_count = DetectAllIntersections();
    if (verbose)
    {
        logoutf("Detected %d intersections", intersection_count);
        logoutf("%u lines after intersection detection", uint32_t(origLines.size()));
    }
    
    std::map<uint32_t, std::set<uint32_t>> collinearOrigLines;

    // sweep all lines
    lines.clear();
    for (auto &line : origLines)
    {
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
                {
                    for (auto &lDup : lines)
                    {
                        if ((lDup.aIdx == aIdx && lDup.bIdx == bIdx) ||
                            (lDup.bIdx == aIdx && lDup.bIdx == bIdx))
                        {
                            logoutf("WARN: Identical line detected! procLine:%u aIdx:P%u bIdx:P%u forOrigLine:#%u", lDup.id, aIdx, bIdx, line.id);
                            collinearOrigLines[lDup.origLine].insert(line.id);
                            foundDup = true;
                        }
                    }
                }
                
                if (foundDup)
                {
                    logoutf("WARN: line #%u has duplicate procLines! aIdx:P%u bIdx:P%u", line.id, aIdx, bIdx);
                }
                else if (!pointsDiffer(intersectionPoints[aIdx], intersectionPoints[bIdx]))
                {
                    logoutf("WARN: line #%u has zero length! aIdx:P%u bIdx:P%u", line.id, aIdx, bIdx);
                }
                else
                {
                    lines.push_back(newLine(i - 1, i, line));
                }
            }
        }
    }
    
    // set collinearity
    // TODO: also check collinearity using area==0, but that's not a real life scenario
    for (auto &kv : collinearOrigLines)
    {
        PolyLine *origLine = nullptr;
        for (auto &l : origLines)
        {
            if (l.id == kv.first)
            {
                origLine = &l;
                break;
            }
        }
        if (!origLine)
        {
            assert(false);
        }
        else
        {
            for (auto &lid : kv.second)
            {
                logoutf("line #%u is collinear with #%u", kv.first, lid);
            }
            for (auto &l : lines)
            {
                if (std::find(kv.second.begin(), kv.second.end(), l.origLine) != kv.second.end())
                {
                    logoutf("set origLine #%u to procLine {%u aIdx:P%u bIdx:P%u}", kv.first, l.id, l.aIdx, l.bIdx);
                    l.origLine = kv.first;
                }
            }
        }
    }
    
    if (verbose)
        logoutf("nOrigLines:%u nLines:%u", uint32_t(origLines.size()), uint32_t(lines.size()));

    return true;
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

/***
* @descr removes line overlappings
* @note must be called before applying Bentley-Ottmann algorithm
*/
void PolyDetector::RemoveOverlappings()
{
#if 0
    size_t i, j, count = origLines.size();
    
#if 1
    for (auto &l : origLines)
    {
        l.calcCenter();
    }
#endif
    
    PolyLine line;

    // lets find overlapping lines
    for (i = 0; i < count; i++)
    {
        auto &line_i = origLines[i];

        for(j = i + 1; j < count; j++)
        {
            auto &line_j = origLines[j];
            if (::Overlap(line_i, line_j))
            {
                int ret = PolyLine::SimplifiedLine(line_i, line_j, line);

                if (ret == 1)
                {
                    // must remove line_j
                    origLines.erase(origLines.begin() + j);
                    j--;
                    count--;
                }
                else
                {
                    if (ret != 2)
                    {
                        // must remove both line_i and line_j and add a new one
                        origLines.erase(origLines.begin() + j);
                        origLines.push_back(line);
                    }

                    // must remove line_i
                    origLines.erase(origLines.begin() + i);

                    // update counters
                    i--;
                    count --;

                    // skip inner loop an go to next step of outer loop
                    break;
                }
            }

        }
    }
#endif
}

uint32_t PolyDetector::DetectAllIntersections()
{
    uint32_t ret = 0;
    size_t counter = origLines.size();
    PointType intersection;
    
    intersectionPoints.clear();
    collinearLineMap.clear();
    
    // sort to always have the same results
    std::sort(origLines.begin(), origLines.end(), PolyLine::bCompareLineOrder);
    
    uint32_t n = 0;
    for (auto &l : origLines)
    {
        l.id = n++;
        l.intersections.clear();
        l.intersectedLines.clear();
        
        l.calcCenter();
    }
    
    // intersected lines: remove lines with only one intersection
    for (uint32_t i = 0; i < counter; i++)
    {
        auto &l1 = origLines[i];
        for (uint32_t j = i + 1; j < counter; j++)
        {
            auto &l2 = origLines[j];
            if (l1.PolyIntersects(l2))
            {
                if (verbose > 1)
                    logoutf("line #%u intersects #%u", l1.id, l2.id);
                
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
    
    if (verbose)
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
                            logoutf("WARN: origLine1 #%u intersects origLine2 #%u exactly in the same point P%u as:", l1.id, l2.id, pi);
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
                            dupPoint = true;
                            intersectionIdx = pi;
                            break;
                        }
                    }

                    if (!dupPoint)
                    {
                        intersectionPoints.push_back(intersection);
                    }
                    
                    if (verbose > 1)
                        logoutf("origLine1 #%u intersects origLine2 #%u in P%u", l1.id, l2.id, intersectionIdx);
                    
                    for (auto l : {&l1, &l2})
                    {
                        if (std::find(l->intersections.begin(), l->intersections.end(), intersectionIdx) != l->intersections.end())
                        {
                            logoutf("WARN: duplicate intersecionIdx P%u in line #%u", intersectionIdx, l->id);
                        }
                        else
                        {
                            l->intersections.push_back(intersectionIdx);
                        }
                    }
                    
                    took.insert(std::make_pair(l1.id, l2.id));
                    took.insert(std::make_pair(l2.id, l1.id));
                    
                    ret++;
                }
            }
        }
    }
    if (verbose)
        logoutf("intersectionPoints:%u", uint32_t(intersectionPoints.size()));
    
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
    str += "]";
    
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
    origLines.push_back(line);
}

void PolyPol::addLine(const PolyLine &l)
{
    p.push_back(l.a);
    p.push_back(l.b);
}

PolyLine *PolyDetector::findLine(uint32_t id, bool useIgnore)
{
    for (auto &l : lines)
    {
        if (l.id == id)
        {
            if (useIgnore && l.ignore)
            {
                //logoutf("line id %u is ignored!", id);
                return nullptr;
            }
            return &l;
        }
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
                
#if 0
                if (l->took > 0 || colStep++ > 1)
                //if (l->took > 0)
                //if (colStep++ > 1)
                {
                    return false;
                }
#else
                return false;
#endif
            }
        }
        
        PolyLine ls(l->center, l1->center);
        for (auto &id2 : idx)
        {
            if (id2 != id && id2 != id1
                //&& !pd.CollinearIdx(id2, id) && !pd.CollinearIdx(id2, id1)
            )
            {
                auto l2 = pd.findLine(id2);
                
                if (l2 && ls.PolyIntersects(*l2))
                //if (l2 && !Collinear(ls, *l2) && ls.PolyIntersects(*l2))
                {
                    if (pd.verbose > 2)
                        logoutf("line %u can't be added to cycle (coll test)! id2:%u intersects with middle(%u, %u)", id, id2, l->id, l1->id);
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
    
    if (cycle.canBeClosed(id))
    {
        cycle.isClosed = true;
        
        if (verbose > 2)
            cycle.print("[CLOSED] ");
        
        if (!similarCycle(_cycles, cycle))
        {
            _cycles.push_back(cycle);
            if (verbose > 2)
                cycle.print("[ACCEPTED] ");
            //logoutf("nCycles:%u", uint32_t(_cycles.size()));
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
    
    if (CycleProcessed(cycle.idx))
    {
        if (verbose > 2)
            cycle.print("[PROCESSED!] ");
        return true;
    }
    processed.insert(cycle.idx);
    
    //if (polyIndices primaryId == polyIndices[0])

    for (auto &nid : _neighbors[id])
    {
        if (_neighbors[nid].size() < 2) continue;
        if (cycle.canBeClosed(nid) || !cycle.contains(nid))
        {
            if (verbose > 2)
                cycle.print((std::string("[NEIGN:") + std::to_string(nid) + std::string("]")).c_str());
            BuildCycle(nid, cycle);
        }
    }
    
    return true;
}

bool PolyDetector::FindPolys()
{
    //logoutf("FindPolys");
    
    processed.clear();
    pointToLines.clear();
    collinearLineMap.clear();
    
    // assign line ids
    uint32_t n = 0;
    for (auto &l : lines)
    {
        //if (l.ignore) continue;
        //logoutf("assign id %u to #%u", n, l.id);
        
        if (dissolveCount == 0) // only on first step
        {
            //logoutf("rename line #%u to %u", l.id, n);
            l.id = n++;
        }
        
        if (!l.ignore)
        {
            pointToLines[l.aIdx].push_back(l.id);
            pointToLines[l.bIdx].push_back(l.id);
            
            //logoutf("line %u a:P%u b:P%u len:%f", l.id, l.aIdx, l.bIdx, l.a.dist2(l.b));
        }
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
                std::sort(neigh.begin(), neigh.end(), [this, &l1](const uint32_t &lid1, const uint32_t &lid2) {
                    auto dl1 = lines[lid1].center.squaredist(l1.center);
                    auto dl2 = lines[lid2].center.squaredist(l1.center);
                    return dl1 < dl2;
                });
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
    
    for (auto &kv : _neighbors) // point by point
    {
        if (_neighbors[kv.first].size() < 2) continue;
        
        if (verbose > 2)
            logoutf("----------- BEGIN L:%u {", kv.first);
        
        PolyCycle cycle;
        cycle.startIdx = kv.first;
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
                
                lPtr->took++;
                
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
    
    if (verbose)
        logoutf("ignore line %u and line %u", l1.id, l2.id);
    l1.setIgnore(*this, "rmCollinear.l1");
    l2.setIgnore(*this, "rmCollinear.l2");

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
    pointToLines[nl.aIdx].push_back(nl.id);
    pointToLines[nl.bIdx].push_back(nl.id);
    
    // build neighbors
    for (auto &id : {l1.id, l2.id})
    {
        for (auto &n : _neighbors[id])
        {
            auto l = findLine(n);
            if (l)
            {
                _neighbors[nl.id].push_back(n);
                _neighbors[n].push_back(nl.id);
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
    
    //return false; // test
    
    lines.push_back(nl);
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
        /*if (nValid > 0)
        {
            logoutf("CollinearIdx[P%u](%s, %s)) = %d. nValid:%u", id, l.toString(*this).c_str(), l1->toString(*this).c_str(), CollinearIdx(l, *l1), nValid);
        }*/
        if (nValid == 1 && CollinearIdx(l, *l1))
        {
            if (verbose)
                logoutf("line %s -> point %s has 2 intersections. otherLine:%u", l.toString(*this).c_str(), id==l.aIdx?"A":"B", l1->id);
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
    
    if (verbose)
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
            assert(false);
            l.took = 2;
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
    for (auto &kv : collinearLineMap)
    {
        if (kv.first == l1)
        {
            return std::find(kv.second.begin(), kv.second.end(), l2) != kv.second.end();
        }
    }
    return false;
}

bool PolyDetector::CollinearIdx(const PolyLine &l1, const PolyLine &l2)
{
    return CollinearIdx(l1.id, l2.id);
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
        std::advance(it2, 1);
        
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
                    
                    if (l3 && ls.PolyIntersects(*l3))
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
    processed.clear();
    origLines.clear();
    lines.clear();
    polys.clear();
    _neighbors.clear();
    collinearLineMap.clear();
    intersectionPoints.clear();
    pointToLines.clear();
    dissolveCount = 0;
}
