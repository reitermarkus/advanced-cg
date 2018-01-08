#include "TriangleMeshLoader.h"

template <typename T>
vector<T> TriangleMeshLoader::readLine(ifstream& file, int maxLine) {
  vector<T> elements;
  string line;
  int lineCount = 0;

  while (getline(file, line)) {
    istringstream iss(line);
    float a, b, c;

    if(lineCount++ == maxLine) {
      break;
    }

    if (!(iss >> a >> b >> c)) {
      elements.push_back(T(0, 0, 0));
    } else {
      elements.push_back(T(a, b, c));
    }
  }

  return elements;
}

optional<vector<Triangle>> TriangleMeshLoader::loadTriangleMesh(string filePath) {
  ifstream file(filePath);
  vector<Triangle> tris;

  if (file) {
    while(file.peek() != EOF) {
      const vector<Vector>& vectors = readLine<Vector>(file, 3);
      const vector<Color>& colors = readLine<Color>(file, 2);

      for (auto& vec : vectors)
        cout << "vec: " << vec << endl;

      for (auto& color : colors)
        cout << "color: " << color << endl;

      cout << endl;

      if(vectors.size() == 3 && colors.size() == 2) {
        tris.push_back(Triangle(vectors[0], vectors[1], vectors[2], colors[0], colors[1], DIFF));
      } else {
        return nullopt;
      }
    }
  } else {
    cerr << "cannot read file" << endl;
  }

  return optional<vector<Triangle>>{tris};
}
