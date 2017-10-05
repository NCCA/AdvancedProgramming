#include "Util.h"

double clamp(double x)
{
  return x < 0 ? 0 : x > 1 ? 1 : x;
}

int toInt(double x)
{
  return int(pow(clamp(x), 1 / 2.2) * 255 + .5);
}

float squareRoot(float x)
{
  unsigned int i = *(unsigned int*) &x;

  // adjust bias
  i  += 127 << 23;
  // approximation of square root
  i >>= 1;

  return *(float*) &i;
}
