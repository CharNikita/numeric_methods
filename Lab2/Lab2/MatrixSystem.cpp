#include "MatrixSystem.h"

#include <fstream>
#include <iostream>
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
   fs.precision(17);
   fs.open(path);

   fs >> n >> m >> k >> max_iter >> block_size;
   fs >> eps;

   di.resize(n);
   for (size_t i = 0; i < di.size(); ++i)
      fs >> di[i];

   al.resize(3);

   al[0].resize(n - 1);
   for (size_t i = 0; i < al[0].size(); ++i)
      fs >> al[0][i];

   al[1].resize(n - m - 2);
   for (size_t i = 0; i < al[1].size(); ++i)
      fs >> al[1][i];

   al[2].resize(n - m - k - 3);
   for (size_t i = 0; i < al[2].size(); ++i)
      fs >> al[2][i];

   au.resize(3);

   au[0].resize(n - 1);
   for (size_t i = 0; i < au[0].size(); ++i)
      fs >> au[0][i];

   au[1].resize(n - m - 2);
   for (size_t i = 0; i < au[1].size(); ++i)
      fs >> au[1][i];

   au[2].resize(n - m - k - 3);
   for (size_t i = 0; i < au[2].size(); ++i)
      fs >> au[2][i];

   b.resize(n);
   for (size_t i = 0; i < b.size(); ++i)
      fs >> b[i];

   x.resize(n);

   normb = norm(b);
   fs.close();
}

// вычисление нормы в эвклидовом порстранстве
template <typename T>
T MatrixSystem<T>::norm(std::vector<T> x)
{
   T result = 0;
   for (int i = 0; i < n; i++)
      result += pow(x[i], 2);

   result = sqrt(result);

   return result;
}

//умножение вектора на вектор
template <typename T>
T MatrixSystem<T>::multVV(int flag, int i, std::vector<T> x0)
{
   T result = 0;
   if (flag == 1 || flag == 3) //умножение нижнего треугольника матрицы
   {
      if (i > 0)
      {
         result += al[0][i - 1] * x0[i - 1];
         if (i > m + 1)
         {
            result += al[1][i - m - 2] * x0[i - m - 2];
            if (i > m + k + 2)
               result += al[2][i - m - k - 3] * x0[i - m - k - 3];
         }
      }
   }

   if (flag == 2 || flag == 3) // умножение верхнего треугольника и диагонали
   {
      result += di[i] * x0[i];
      if (i < n - 1)
      {
         result += au[0][i] * x0[i + 1];
         if (i < n - m - 2)
         {
            result += au[1][i] * x0[i + m + 2];
            if (i < n - m - k - 3)
               result += au[2][i] * x0[i + m + k + 3];
         }
      }
   }

   return result;
}

// метод якоби
template <typename T>
T MatrixSystem<T>::jacobi(T w, std::vector<T>& x1, T& loss)
{
   T sum = 0;
   T buf = 0;
   loss = 0;

   for (int i = 0; i < n; i++)
   {
      sum = multVV(3, i, x);
      buf = b[i] - sum;
      x1[i] = x[i] + w * buf / di[i];
      loss += pow(buf, 2);
   }
   loss = sqrt(loss) / normb;

   return loss;
}

template <typename T>
T MatrixSystem<T>::gaussZeidel(T w, std::vector<T>& x1, T& loss)
{
   T sum = 0;
   T buf = 0;
   loss = 0;

   for (int i = 0; i < n; i++)
   {
      sum = multVV(2, i, x) + multVV(1, i, x1);
      buf = b[i] - sum;
      x1[i] = x[i] + (w / di[i]) * buf;
      loss += pow(buf, 2);
   }
   loss = sqrt(loss) / normb;

   return loss;
}

// подсчет числа обусловленностей для метода якоби и гаусса-зейделя
template <typename T>
T MatrixSystem<T>::num_obusl(std::vector<T> x, T loss, T normxstar)
{
   T obusl;
   for (size_t i = 0; i < n; i++)
   {
      x[i] -= i + 1;
   }
   obusl = norm(x);
   obusl = (obusl / normxstar) / loss;
   return obusl;

}

// метод якоби (flag = false) и гаусса-зейделя (flag = true)
template <typename T>
void MatrixSystem<T>::iteration(std::string path, bool flag)
{
   T obusl = 0;
   T w = 0;
   T loss = 0;
   T normxstar = 0;
   std::ofstream fs;
   //std::locale mylocale("");
   fs.open("output.txt");
   fs.imbue(std::locale("Russian"));
   fs.precision(17);

   size_t  t;
   for (size_t i = 1; i <= n; i++)
      normxstar += pow(i, 2);

   normxstar = sqrt(normxstar);
   std::vector<T> buf(n);
   for (int i = 171; i <= 175; i += 1)
   {
      w = i / 100.0;
      int exit = 0;
      std::cout << w << std::endl;
      for (t = 0; t < max_iter && exit != 1; t++)
      {
         if (flag == false)
         {
            jacobi(w, buf, loss);
            x.swap(buf);
         }

         if (flag == true)
         {
            gaussZeidel(w, buf, loss);
            x.swap(buf);
         }

         if (loss < eps)
         {
            exit = 1;
            std::cout << "End iter:" << t << std::endl;
            std::cout << "EXIT" << std::endl;
         }

         obusl = num_obusl(x, loss, normxstar);
         if (t % 1000 == 0)
         {
            std::cout << "iter: " << t + 1 << std::endl;
            std::cout << "loss: " << loss << std::endl;
            std::cout << "obusl: " << obusl << std::endl;
         }
      }

      fs << w << "\t" << x[0] << "\t" << t + 1 << "\t" << loss << "\t" << obusl << std::endl;
      for (size_t j = 1; j < n; j++)
         fs << "\t" << x[j] << std::endl;
      fs << std::endl;

      for (size_t j = 0; j < n; j++)
      {
         x[j] = 0;
         buf[j] = 0;
      }
      loss = 0;
   }

   fs.close();
}


template <typename T>
void MatrixSystem<T>::lu()
{
   //int block_count = n / block_size;
   for (int i = 0; i < n / block_size; i++)
   {
      int k = block_size - 1 + i * block_size;
      for (int j = i * block_size; j < k; j++)
      {
         au[0][j] = au[0][j] / di[j];
         di[j + 1] -= al[0][j] * au[0][j];
      }
   }
}

template <typename T>
void MatrixSystem<T>::output()
{
   std::ofstream fs;
   fs.open("output.txt");

   for (size_t i = 0; i < n; i++)
      fs << std::fixed << x[i] << "\n";

   fs.close();
}
