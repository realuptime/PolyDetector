
#pragma once

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
};

