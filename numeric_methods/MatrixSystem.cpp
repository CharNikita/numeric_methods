#include "MatrixSystem.h"

#include <fstream>
#include <string>

template <typename T>
MatrixSystem<T>::MatrixSystem(std::string path)
{
   readFromFile(path);
}

template<typename T>
MatrixSystem<T>::~MatrixSystem()
{
   delete di;
   delete ia;
   delete al;
   delete au;
   delete b;
}


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

   al.resize(size + 1);
   for (size_t i = 0; i < size + 1; ++i)
      fs >> al[i];
   
   au.resize(size + 1);
   for (size_t i = 0; i < size + 1; ++i)
      fs >> au[i];
   
   b.resize(size);
   for (size_t i = 0; i < size; ++i)
      fs >> b[i];

   fs.close();
}
