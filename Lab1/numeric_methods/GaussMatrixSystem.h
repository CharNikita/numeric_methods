#pragma once

#include <vector>
#include <string>

template <typename T>
class GaussMatrixSystem
{
   std::vector<std::vector<T>> M;
   std::vector<T> b;
   size_t size;

   void readFromFile(std::string path);
   GaussMatrixSystem();

public:
   size_t getSize() { return size; }

   GaussMatrixSystem(std::string path);

   ~GaussMatrixSystem();

   void gauss();
   void output();
};