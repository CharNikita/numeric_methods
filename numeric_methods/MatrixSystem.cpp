#include "MatrixSystem.h"

#include <fstream>
#include <string>

template <typename T>
MatrixSystem<T>::MatrixSystem(std::string path)
{
   readFromFile(path);
}

template <typename T>
void MatrixSystem<T>::readFromFile(std::string path)
{
   std::fstream fs;
   fs.open(path);

   fs >> size;
   for (size_t i = 0; i < size; ++i)
      fs >> di[i];

   for (size_t i = 0; i < size; ++i)
      fs >> ia[i];

   for (size_t i = 0; i < size; ++i)
      fs >> al[i];

   for (size_t i = 0; i < size; ++i)
      fs >> au[i];

   for (size_t i = 0; i < size; ++i)
      fs >> b[i];

   fs.close();
}
