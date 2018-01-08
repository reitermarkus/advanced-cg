#ifndef __TRIANGLE_MESH_LOADER_H__
#define __TRIANGLE_MESH_LOADER_H__

#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <optional>

#include "Triangle.h"
#include "Vector.h"

using namespace std;

class TriangleMeshLoader {
  private:
    template <typename T>
    static vector<T> readLine(ifstream& file, int maxLine);

  public:
    static optional<vector<Triangle>> loadTriangleMesh(string filePath);
};

#endif
