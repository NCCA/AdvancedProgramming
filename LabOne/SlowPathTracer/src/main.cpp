// base on smallpt, a Path Tracer by Kevin Beason, 2008s

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <array>
#include <algorithm>
#include <iostream>
#include <random>
#include <memory>
#include <chrono>
#include <limits>
#include "Vec.h"
#include "Ray.h"
#include "Sphere.h"
#include "Material.h"
#include "Util.h"
#include "Framebuffer.h"
#include "Image.h"

static std:: default_random_engine generator;
static std::uniform_real_distribution<double> distr(0.0,1.0);
double uniformRand()
{
    return distr(generator);
}



static std::array<Sphere,9> spheres = {{
  //Scene: radius, position, emission, color, material
  Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.95, .25, .25), DIFF), //Left
  Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .95), DIFF), //Rght
  Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF),   //Back
  Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFF),         //Frnt
  Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF),   //Botm
  Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF), //Top
  Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1) * .999, SPEC),  //Mirr
  Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1) * .999, REFR),  //Glas
  Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(12, 12, 12), Vec(), DIFF) //Lite
}};


bool intersect(const Ray& r, double& t, size_t& id)
{
  double d;
  double inf=t=std::numeric_limits<double>::max();

  for (size_t i = spheres.size(); i--; )
  {
    if ((d = spheres[i].intersect(r)) && d < t)
    {
      t  = d;
      id = i;
    }
  }
  return t < inf;
}

Vec radiance(const Ray& r, int depth)
{
  double t;                         // distance to intersection
  size_t id = 0;                       // id of intersected object
  if (!intersect(r, t, id))
  {
    return Vec(); // if miss, return black
  }
  const Sphere& obj = spheres[id];  // the hit object
  Vec x = r.o + r.d * t;
  Vec n = (x - obj.p).norm();
  Vec nl = n.dot(r.d) < 0 ? n : n * -1;
  Vec f = obj.c;
  double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;  // max refl
  if (++depth > 5)
  {
    if (uniformRand() < p)
    {
      f = f * (1 / p);
    }
    else
    {
      return obj.e; //R.R.
    }
  }
  if (obj.refl == DIFF)
  { // Ideal DIFFUSE reflection
    double r1 = 2 * M_PI * uniformRand(), r2 = uniformRand(), r2s = sqrt(r2);
    Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(), v = w % u;
    Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
    return obj.e + f.mult(radiance(Ray(x, d), depth));
  }
  else if (obj.refl == SPEC)
  { // Ideal SPECULAR reflection
    return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth));
  }
  Ray reflRay(x, r.d - n * 2 * n.dot(r.d)); // Ideal dielectric REFRACTION
  bool into = n.dot(nl) > 0;             // Ray from outside going in?
  double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;
  if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Total internal reflection
  {
    return obj.e + f.mult(radiance(reflRay, depth));
  }
  Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
  double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.dot(n));
  double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
  return obj.e + f.mult(depth > 2 ? (uniformRand() < P ? // Russian roulette
                                     radiance(reflRay, depth) * RP
             : radiance(Ray(x, tdir), depth) * TP)
            : radiance(reflRay, depth) * Re + radiance(Ray(x, tdir), depth) * Tr);
}





int main(int argc, char* argv[])
{


  size_t w =400;
  size_t h = 400;
  size_t samps = 4;
  int opt;
  std::string fname("image.ppm");
  // quick and dirty command line check no error checking!
  while ((opt = getopt(argc,argv,"w:h:s:f:")) != EOF)
  {
    switch(opt)
    {
      case 'w' : w=atoi(optarg); break;
      case 'h' : h=atoi(optarg); break;
      case 's' : samps=atoi(optarg); break;
      case 'f' : fname=optarg; break;
      case '?' :
        std::cerr<<"usage PathTracer -w [width] -h [height] -s [samples] \n";
        exit(EXIT_SUCCESS);
      break;
    }
  }

#ifdef USEFRAMEBUFFER

  std::unique_ptr<frm::Framebuffer> framebuffer( new frm::Framebuffer());
  framebuffer->init(w, h, NULL);

  Image image(w,h);
  image.setBackground(255,255,255);

  framebuffer->bind();
  framebuffer->poll();
  framebuffer->title("smallpt");
#endif
  //size_t cPixel=0;

  Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm());        // cam pos, dir
  Vec cx = Vec(w * .5135 / h);
  Vec cy = (cx % cam.d).norm() * .5135;
  Vec r;
  std::unique_ptr<Vec[]> c(new Vec[w * h]);
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();

  for (size_t y = 0; y < h; y++)
  {                    // Loop over image rows
    #ifdef USEFRAMEBUFFER
        char msg[50];

        sprintf(msg, "Rendering (%ld spp) %5.2f%%", samps * 4, 100. * y / (h - 1));
        framebuffer->title(std::string(msg));
    #else
      fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
    #endif
    for (unsigned short x = 0; x < w; x++)
    {// Loop cols
      for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++) // 2x2 subpixel rows
      {
        for (int sx = 0; sx < 2; sx++, r = Vec())
        { // 2x2 subpixel cols
          for (size_t s = 0; s < samps; s++)
          {
            double r1 = 2 * uniformRand(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
            double r2 = 2 * uniformRand(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
            Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
            cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
            r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0) * (1. / samps);
          } // Camera rays are pushed ^^^^^ forward to start in interior
          c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
        }
        #ifdef USEFRAMEBUFFER
          image.setPixel(x,h-y,toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
        #endif
      }

    }
    #ifdef USEFRAMEBUFFER
      framebuffer->image(image.get(), w, h);
      framebuffer->poll();
      if(framebuffer->shouldClose())
      {
        exit(EXIT_FAILURE);
      }
      framebuffer->draw();
    #endif
    } // end for rows

  end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  std::cout << "finished computation at " << std::ctime(&end_time)
                << "elapsed time: " << elapsed_seconds.count() << "s\n";

#ifdef USEFRAMEBUFFER
  while(!framebuffer->shouldClose())
  {
    framebuffer->poll();
    framebuffer->draw();
    sleep(1);
  }
#endif
  FILE *f = fopen(fname.c_str(), "w");         // Write image to PPM file.
  fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
  for (int i=0; i<w*h; i++)
  {
    fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
  }
  fclose(f);

}
