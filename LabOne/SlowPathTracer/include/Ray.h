#ifndef RAY_H_
#define RAY_H_
#include "Vec.h"
struct Ray
{
  Vec o, d;
  Ray(Vec o_, Vec d_)
    : o(o_), d(d_) {
  }
};

#endif
