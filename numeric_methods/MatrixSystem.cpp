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

      // k - index of current element in 'al' and 'au'
      for (int k = ia[i] - 1; k < ia[i + 1] - 1; ++k)
      {
         // sums for j-th row/column
         T sum_al = 0;
         T sum_au = 0;

         // kj - index of first element in a row for L
         //                      and in a column for U
         int kj = ia[i] - 1;

         // ki - index of first element in a column for L
         //                            and in a row for U
         int ki = ia[k - kj] - 1 + i - count;

         // index of first diagonal element
         int kd = i - count;

         // difference between numbers of elements
         //    for L: num row's elems - num column's elems
         //    for U: num column's elems - num row's elems
         const int diff = (k + 1 - ia[i]) - 
            (ia[k - kj + i - count + 1] - ia[k - kj + i - count]);

         // if 'diff' is positive then increase indexes 'kj' and 'kd'
         // else increase index 'ki'
         if (diff > 0)
         {
            kj += diff;
            kd += diff;
         }
         else
            ki += -diff;

         for (kj; kj < k; ++kj)
         {
            sum_al += al[kj] * di[kd] * au[ki];
            sum_au += au[kj] * di[kd] * al[ki];
            ++ki;
            ++kd;
         }

         if (di[kd] == 0) throw std::runtime_error("Matrix cannot be factorized!");

         au[k] = (au[k] - sum_au) / di[kd];
         al[k] = (al[k] - sum_al) / di[kd];
         sum += al[k] * di[kd] * au[k];
      }

      di[i] -= sum;      
   }
}
