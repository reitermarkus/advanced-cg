#ifndef __TRIANGLE_MESH_LOADER_H__
#define __TRIANGLE_MESH_LOADER_H__

#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <map>
#include <optional>

#include "Triangle.h"
#include "Vector.h"

using namespace std;

class TriangleMeshLoader {
  private:
    template <typename T>
    static vector<T> readLine(ifstream& file, int maxLine);
    static void checkReflection(ifstream& file, Refl_t& refl);

  public:
    static optional<vector<Triangle>> loadTriangleMesh(string filePath);
};

struct ReflMap : public map<string, Refl_t> {
  ReflMap() {
    this->operator[]("DIFF") = DIFF;
    this->operator[]("SPEC") = SPEC;
    this->operator[]("REFR") = REFR;
  };
};

#endif
