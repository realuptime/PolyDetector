
#pragma once

struct vec
{
    float x, y, z;

    vec() {}

    vec(float a, float b, float c) : x(a), y(b), z(c) {}
    vec(float a) : x(a), y(a), z(a) {}

    vec &sub(const vec &v)
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }
    vec &add(const vec &v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    vec &mul(const float &v)
    {
        x *= v;
        y *= v;
        z *= v;
        return *this;
    }
    vec &div(const float &v)
    {
        x /= v;
        y /= v;
        z /= v;
        return *this;
    }

    float squaredlen() const
    {
        return x*x + y*y + z*z;
    }

    float squaredist(const vec &v) const
    {
        return vec(*this).sub(v).squaredlen();
    }
};

