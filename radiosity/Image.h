#ifndef __IMAGE_H__
#define __IMAGE_H__

#include <cstring>
#include <fstream>
#include <iostream>

#include "Vector.h"

struct Image {
  size_t width, height;
  Color *pixels;

  Image(size_t _w, size_t _h);

  Color getColor(size_t x, size_t y);
  void setColor(size_t x, size_t y, const Color &c);
  void addColor(size_t x, size_t y, const Color &c);
  void save(const string &filename);
};

#endif //__IMAGE_H__
