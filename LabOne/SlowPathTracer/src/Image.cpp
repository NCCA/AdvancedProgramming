#include "Image.h"
constexpr int c_pixelDepth=3;
Image::Image(size_t _width, size_t _height)
{
  m_data.reset(new unsigned char[_width*_height*c_pixelDepth]);
  m_width=_width;
  m_height=_height;
}

void Image::setPixel(size_t _x, size_t _y, unsigned char _r,  unsigned char _g,unsigned char _b,unsigned char _a)
{
  // ignore out of range pixels
  if( (_x <=0 || _x>=m_width) || (_y<=0 || _y>=m_height) )
    return;
  size_t offset=(_y*c_pixelDepth * m_width) + _x*c_pixelDepth;
  m_data[offset]=_r;
  m_data[offset+1]=_g;
  m_data[offset+2]=_b;

}



void Image::setBackground(unsigned char _r, unsigned char _g, unsigned char _b)
{
  for(size_t y=0; y<m_height; ++y)
  {
    for(size_t x=0; x<m_width; ++x)
   {
      setPixel(x,y,_r,_g,_b);
    }
  }
}
