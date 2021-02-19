
#pragma once

#include <limits>
#include <cmath>

struct vec
{
    float x, y;

    vec() {}

    vec(float a, float b) : x(a), y(b) {}
    vec(float a) : x(a), y(a) {}
    vec(float a, float b, float c) : x(a), y(b) {}

    vec &sub(const vec &v)
    {
        x -= v.x;
        y -= v.y;
        return *this;
    }
    vec &add(const vec &v)
    {
        x += v.x;
        y += v.y;
        return *this;
    }
    vec &mul(const float &v)
    {
        x *= v;
        y *= v;
        return *this;
    }
    vec &div(const float &v)
    {
        x /= v;
        y /= v;
        return *this;
    }

    float squaredlen() const
    {
        return x*x + y*y;
    }

    float squaredist(const vec &v) const
    {
        return vec(*this).sub(v).squaredlen();
    }
	
	    static void line(const vec &p, const vec &q, float &a, float &b, float &c)
    {
        // Line AB represented as a*x + b*y = c
        a = p.y - q.y;
        b = q.x - p.x;
        c = (p.x - q.x)*p.y + (q.y - p.y)*p.x;
    }
    
    static float lineDist(float a, float b, float c, const vec &p)
    {
        if (a + b == 0.0f) return std::numeric_limits<float>::max();;
        return fabs(a * p.x + b * p.y + c) / sqrtf(a * a + b * b);
    }
    
    static float lineDist(const vec &la, const vec &lb, const vec &p)
    {
        float a, b, c;
        line(la, lb, a, b, c);
        return lineDist(a, b, c, p);
    }
};

