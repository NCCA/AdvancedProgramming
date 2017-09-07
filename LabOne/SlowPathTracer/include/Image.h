#ifndef IMAGE_H_
#define IMAGE_H_
#if defined (__linux__) || defined (WIN32)
  #include <GL/gl.h>
  #include <GL/glu.h>
#endif
#ifdef __APPLE__
  #include <OpenGL/gl.h>
  #include <OpenGL/glu.h>
#endif
#include <memory>
#include <string>

class Image
{
  public :

    explicit Image(size_t _width, size_t _height);
    // make sure we can't copy this class as it could be expensive.
    Image(const Image &)=delete;
    Image & operator=(const Image &)=delete;
    void setPixel(size_t _x, size_t _y, unsigned char _r,
                                        unsigned char _g,
                                        unsigned char _b, unsigned char _a=1.0f);
    void setBackground(unsigned char _r,
                       unsigned char _g,
                       unsigned char _b);
    unsigned char* get() {return m_data.get();}
  private :

   std::unique_ptr< unsigned char []> m_data;
   size_t m_width=0;
   size_t m_height=0;
};


#endif
