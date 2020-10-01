#include "GaussMatrixSystem.h"

#include <fstream>
#include <string>

template <typename T>
GaussMatrixSystem<T>::GaussMatrixSystem(std::string path)
{
   readFromFile(path);
}

template <typename T>
GaussMatrixSystem<T>::~GaussMatrixSystem() {}

template <typename T>
void GaussMatrixSystem<T>::readFromFile(std::string path)
{
   std::ifstream fs;
   fs.open(path);

   fs >> size;

   M.resize(size, std::vector<T>(size));
   for (size_t i = 0; i < size; ++i)
      for (size_t j = 0; j < size; ++j)
         fs >> M[i][j];

   b.resize(size);
   for (size_t i = 0; i < size; ++i)
      fs >> b[i];

   fs.close();
}

template <typename T>
void GaussMatrixSystem<T>::gauss()
{
   for (int i = 1; i < size; ++i)
      for (int j = i; j < size; ++j)
      {
         if (M[i - 1][i - 1] == 0) throw std::runtime_error("Result cannot be calculated!");

         if (M[j][i - 1] != 0)
         {

            T mult = -M[j][i - 1] / M[i - 1][i - 1];

            for (int k = i - 1; k < size; ++k)
            {
               M[j][k] += M[i - 1][k] * mult;
            }

            b[j] += b[i - 1] * mult;
         }
      }

   for (int i = size - 1; i >= 0; --i)
   {
      T sum = 0;
      for (int j = size - 1; j > i; --j)
      {
         sum += M[i][j] * b[j];
      }

      b[i] = (b[i] - sum) / M[i][i];
   }
}

template <typename T>
void GaussMatrixSystem<T>::output()
{
   std::ofstream fs;
   fs.open("output.txt");

   for (size_t i = 0; i < size; i++)
      fs << std::fixed << b[i] << "\n";

   fs.close();
}
