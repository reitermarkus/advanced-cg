#include "TriangleMeshLoader.h"

template <typename T>
vector<T> TriangleMeshLoader::readLine(ifstream& file, int maxLine) {
  vector<T> elements;
  string line;

  for (int lineCount = 0; lineCount < maxLine; lineCount++) {
    getline(file, line);
    istringstream iss(line);
    float a, b, c;

    if (!(iss >> a >> b >> c)) {
      elements.push_back(T(0, 0, 0));
    } else {
      elements.push_back(T(a, b, c));
    }
  }

  return elements;
}

void TriangleMeshLoader::checkReflection(ifstream& file, Refl_t& refl) {
  string reflectionString;

  auto rightTrim = [] (string &s) {
    s.erase(find_if(s.rbegin(), s.rend(), [] (int ch) {
      return !isspace(ch);
    }).base(), s.end());

    return s;
  };

  if(getline(file, reflectionString)) {
    istringstream iss(reflectionString);
    float a, b, c;

    if (!(iss >> a >> b >> c)) {
      ReflectionMap reflectionMap;

      if(reflectionMap.find(rightTrim(reflectionString)) != reflectionMap.end()) {
        refl = reflectionMap[reflectionString];
      }
    } else {
      file.clear();
      file.seekg(0, ios::beg);
    }
  }
}

optional<vector<Triangle>> TriangleMeshLoader::loadTriangleMesh(string filePath) {
  ifstream file(filePath);
  vector<Triangle> tris;
  Refl_t reflection = DIFF;

  if (file) {
    checkReflection(file, reflection);

    while(file.peek() != EOF) {
      const vector<Vector>& vectors = readLine<Vector>(file, 3);
      const vector<Color>& colors = readLine<Color>(file, 2);

      if(vectors.size() == 3 && colors.size() == 2) {
        tris.push_back(Triangle(vectors[0], vectors[1], vectors[2], colors[0], colors[1], reflection));
      } else {
        cerr << "invalid Triangle in file" << endl;
        file.close();
        return nullopt;
      }

      file.ignore(numeric_limits<streamsize>::max(), file.widen('\n'));
    }

    file.close();
  } else {
    cerr << "cannot read file" << endl;
    file.close();
    return nullopt;
  }

  return optional<vector<Triangle>>{tris};
}
