#include "Image.h"

#include <algorithm>
#include <cstddef>

using namespace std::literals::string_literals;

Image::Image(size_t _w, size_t _h) : width(_w), height(_h) {
  pixels = new Color[width * height];
}

Color Image::getColor(size_t x, size_t y) {
  size_t image_index = (height - y - 1) * width + x;
  return pixels[image_index];
}

void Image::setColor(size_t x, size_t y, const Color &c) {
  size_t image_index = (height - y - 1) * width + x;
  pixels[image_index] = c;
}

void Image::addColor(size_t x, size_t y, const Color &c) {
  size_t image_index = (height - y - 1) * width + x;
  pixels[image_index] = pixels[image_index] + c;
}

int toInteger(double x) { 
	/* Clamp to [0,1] */
	if (x<0.0) x = 0.0;		
        
	if (x>1.0) x = 1.0;             

	/* Apply gamma correction and convert to integer */
	return int(pow(x,1/2.2)*255+.5); 
}

void Image::save(const string &filename) {
	/* Save image in PPM format */
	FILE *f = fopen(filename.c_str(), "wb");
	fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
	for (int i = 0; i < width * height; i++)
	fprintf(f,"%d %d %d ", toInteger(pixels[i].x), 
							toInteger(pixels[i].y), 
							toInteger(pixels[i].z));
	fclose(f);
}
/*
static std::byte toByte(double x) {
  x = clamp(x, 0.0, 1.0);

  // Apply gamma correction and convert to integer.
  return (byte)int(pow(x, 1.0 / 2.2) * 255.0 + 0.5);
}

// Save image in binary PPM format.
void Image::save(const string &filename) {
  ofstream file;

  if (filename.find(".ppm") != string::npos)
    file.open(filename, ios::out | ios::binary);
  else
    file.open(filename + ".ppm"s, ios::out | ios::binary);

  file << "P6" << endl;
  file << width << " " << height << endl;
  file << 255 << endl;
  for (size_t i = 0; i < width * height; i++) {
    file << (unsigned char)toByte(pixels[i].x);
    file << (unsigned char)toByte(pixels[i].y);
    file << (unsigned char)toByte(pixels[i].z);
  }
  file.close();
}
*/
