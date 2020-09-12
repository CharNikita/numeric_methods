#include "MatrixSystem.h"

#include <fstream>
#include <string>
#include <iostream>

template <typename T>
MatrixSystem<T>::MatrixSystem(std::string path)
{
   readFromFile(path);
}

template<typename T>
MatrixSystem<T>::~MatrixSystem() {}

template <typename T>
void MatrixSystem<T>::readFromFile(std::string path)
{
   std::fstream fs;
   fs.open(path);

   fs >> size;

   di.resize(size);
   for (size_t i = 0; i < size; ++i)
      fs >> di[i];

   ia.resize(size + 1);
   for (size_t i = 0; i < size + 1; ++i)
      fs >> ia[i];

   al.resize(ia[size] - 1);
   for (size_t i = 0; i < ia[size] - 1; ++i)
      fs >> al[i];
   
   au.resize(ia[size] - 1);
   for (size_t i = 0; i < ia[size] - 1; ++i)
      fs >> au[i];
   
   b.resize(size);
   for (size_t i = 0; i < size; ++i)
      fs >> b[i];

   fs.close();
}

template <typename T>
void MatrixSystem<T>::ldu()
{
   for (int i = 1; i < size; ++i)
   {
      // sum for i-th row/column
      T sum = 0;
      // number of items for i-th row/column
      const int count = ia[i + 1] - ia[i];
      // index of current element in 'al' and 'au'
      int k = ia[i] - 1;

      for (int j = i - count; j < i; ++j)
      {
         // sums for j-th row/column
         T sum_al = 0;
         T sum_au = 0;

         int i_row = ia[i + 1] - ia[i];
         int j_row = ia[j + 1] - ia[j];
         int min = j_row < i_row
            ? j_row - (i - count)
            : i_row - (i - count);

         int ki = k - min;
         int kd = k - min - 1 - (i - count);
         for (int kj = k - min - 1; kj < k - min + j - 1 - (i - count); ++kj)
         {
            sum_al += au[kj] * di[kd] * al[ki];
            sum_au += al[kj] * di[kd] * au[ki];
            ++kd;
            ++ki;
         }

         if (di[j] == 0) throw std::runtime_error("Matrix cannot be factorized!");

         au[k] = (au[k] - sum_au) / di[j];
         al[k] = (al[k] - sum_al) / di[j];
         sum += al[k] * di[j] * au[k];
         ++k;
      }

      di[i] -= sum;      
   }
}
