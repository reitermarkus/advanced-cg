#ifndef __IMAGE_H__
#define __IMAGE_H__

#include <cstring>
#include <fstream>
#include <iostream>

#include "Vector.h"

struct Image {
  int width, height;
  Color *pixels;

  Image(int _w, int _h);

  Color getColor(int x, int y);
  void setColor(int x, int y, const Color &c);
  void addColor(int x, int y, const Color &c);
  int toInteger(double x);
  void Save(const string &filename);
};

#endif //__IMAGE_H__
