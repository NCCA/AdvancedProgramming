#ifndef VEC_H_
#define VEC_H_

#include <cmath>

struct Vec
{

  double x=0;
  double y=0;
  double z=0;

  Vec(const Vec& )=default;
  Vec& operator=(const Vec& )=default;
  Vec(Vec&& )=default;
  Vec& operator=(Vec&& )=default;
  ~Vec()=default;
  Vec(double _x = 0, double _y = 0, double _z = 0)
  {
    x = _x;
    y = _y;
    z = _z;
  }

  Vec operator+(const Vec& b) const
  {
    return Vec(x + b.x, y + b.y, z + b.z);
  }
  Vec operator-(const Vec& b) const
  {
    return Vec(x - b.x, y - b.y, z - b.z);
  }
  Vec operator*(double b) const
  {
    return Vec(x * b, y * b, z * b);
  }
  Vec mult(const Vec& b) const
  {
    return Vec(x * b.x, y * b.y, z * b.z);
  }
  Vec& norm()
  {
    return *this = *this * (1 / sqrt(x * x + y * y + z * z));
  }
  double dot(const Vec& b) const
  {
    return x * b.x + y * b.y + z * b.z;
  }                                                                // cross:

  Vec operator%(Vec& b)
  {
    return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
  }
};

#endif
