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
   for (size_t i = 0; i < size; ++i)
   {
      // sum for i-th row/column
      T sum = 0;
      // number of items for i-th row/column
      const int count = ia[i + 1] - ia[i];
      // start position for 'al' and 'au'
      int k = ia[i] - 1;

      for (int j = i - count; j < i; ++j)
      {
         // sums for j-th row/column
         T sum_al = 0;
         T sum_au = 0;
         // current position for 'al' and 'au'
         int numj = ia[j + 1] - ia[j];

         const int max = numj < k - ia[i] ?
            numj - 1 : k - ia[i] - 1;
         
         int fku = ia[j + 1] - 1;
         // current position for 'al' for sum
         int kl = k - 1;
         // current position for 'al'
         int kd = j - 1;
         for (int ku = fku; ku > fku - max; --ku)
         {
            sum_al += au[ku] * di[kd] * al[kl];
            sum_au += al[ku] * di[kd] * au[kl];
            --kl;
            --kd;
         }

         if (di[j] == 0) throw std::runtime_error("ERROR: Could not manipulate with this matrix!");

         au[k] = (au[k] - sum_au) / di[j];
         al[k] = (al[k] - sum_al) / di[j];
         sum += al[k] * di[j] * au[k];
         ++k;
      }

      di[i] -= sum;

      /*if (ia[i + 1] > ia[i])
      {
         sum = 0;
         for (size_t j = ia[i]; j < ia[i + 1]; ++j)
         {
            // sum +=
            T l = al[j];
            T u = au[j];
            l = l - sum;

         }
      }*/
      
   }
}
