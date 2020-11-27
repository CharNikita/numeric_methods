#include "MatrixSystem.h"

#include <fstream>
#include <iostream>
#include <string>
#include <ctime>

template <typename T>
MatrixSystem<T>::MatrixSystem(std::string& path)
{
   readFromFile(path);
}

template <typename T>
void MatrixSystem<T>::readFromFile(std::string& path)
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
   buf.resize(n);

   normb = norm(b);
   fs.close();
}


template <typename T>
std::vector<T> MatrixSystem<T>::get_x()
{
   return x;
}

template <typename T>
void MatrixSystem<T>::set_x(std::vector<T> x)
{
   this->x = x;
}

// вычисление нормы в эвклидовом порстранстве
template <typename T>
T MatrixSystem<T>::norm(std::vector<T>& x)
{
   T result = 0;
   for (int i = 0; i < n; i++)
      result += x[i] * x[i];

   result = sqrt(result);

   return result;
}

//умножение вектора на вектор
template <typename T>
T MatrixSystem<T>::multVV(int flag, int i, std::vector<T>& x0)
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
      x1[i] = x[i] + (w / di[i]) * buf;
      loss += buf * buf;
   }
   loss = sqrt(loss) / normb;

   return loss;
}

// метод гаусса-зейделя
template <typename T>
T MatrixSystem<T>::gauss_zeidel(T w, T& loss)
{
   T sum = 0;
   T buf = 0;
   loss = 0;

   for (int i = 0; i < n; i++)
   {
      sum = multVV(2, i, x) + multVV(1, i, x);
      buf = b[i] - sum;
      x[i] = x[i] + (w / di[i]) * buf;
      loss += buf * buf;
   }
   loss = sqrt(loss) / normb;

   return loss;
}

// подсчет числа обусловленностей для метода якоби и гаусса-зейделя
template <typename T>
T MatrixSystem<T>::num_obusl(std::vector<T> x, T loss, T normxstar)
{
   T obusl;
   for (int i = 0; i < n; i++)
      x[i] -= i + 1;

   obusl = norm(x);
   obusl = obusl / normxstar / loss;

   return obusl;
}

// метод якоби (flag = false) и гаусса-зейделя (flag = true)
template <typename T>
void MatrixSystem<T>::iteration(std::string& path, bool flag)
{
   T obusl = 0;
   T w = 0;
   T loss = 0;
   T normxstar = 0;
   std::ofstream fs;
   fs.open(path);
   fs.imbue(std::locale("Russian"));
   fs.precision(17);

   size_t  t;
   for (size_t i = 1; i <= n; i++)
      normxstar += i * i;

   normxstar = sqrt(normxstar);
   for (int i = 170; i <= 175; i += 1)
   {
      clock_t start = clock();
      w = i / 100.0;
      std::cout << w << std::endl;
      for (t = 0; t < max_iter; t++)
      {         
         if (flag)
            gauss_zeidel(w, loss);
         else
         {
            jacobi(w, buf, loss);
            x.swap(buf);
         }

         if (loss < eps)
         {
            std::cout << "End iter:" << t << std::endl;
            clock_t finish = clock();
            std::cout << "Time: " << finish - start << std::endl;
            std::cout << "EXIT" << std::endl;
            break;
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
void MatrixSystem<T>::lu_factorization()
{
   for (size_t i = 0; i < n / block_size; i++)
   {
      int k = block_size - 1 + i * block_size;
      for (int j = i * block_size; j < k; j++)
      {
         au[0][j] = au[0][j] / di[j];
         di[j + 1] -= al[0][j] * au[0][j];
      }
   }
}

// решение слау прямым и обратным ходом для блока
template <typename T>
void MatrixSystem<T>::lu_solution(std::vector<T>& x0, int block_amt, T w)
{
   int j = (block_amt + 1) * block_size;
   for (size_t i = block_amt * block_size; i < j; i++)
      x0[i] *= w;
   int c = block_amt * block_size;
   x0[c] /= di[c];
   for (size_t i = c + 1; i < j; i++)
   {
      x0[i] -= x0[i - 1] * al[0][i - 1];
      if (di[i] == 0.0)
         throw 2;
      x0[i] /= di[i];
   }
   for (int i = j - 2; i >= j - block_size && i > -1; i--)
   {
      x0[i] -= au[0][i] * x0[i + 1];
   }

}

// метод блочной релаксации: одна итерация
template <typename T>
void MatrixSystem<T>::block_iteration(T w, std::vector<T>& xkp, T& loss_bl)
{
   T sum = 0;
   T buf = 0;
   loss_bl = 0;
   xkp = x;
   for (int i = 0; i < n / block_size; i++)
   {
      multBV(i, xkp); // sum Aki*Xk k = 1 to block_amt 
      int end = (i + 1) * block_size;
      for (int k = i * block_size; k < end; k++)
      {
         xkp[k] = b[k] - xkp[k]; // F - sum Aki*Xk k = 1 to block_amt 
      }
      lu_solution(xkp, i, w); // AiiYi = xkpi 
      for (int k = i * block_size; k < end; k++)
      {
         xkp[k] += (1 - w) * x[k]; // Xk+1 = Yk + (1 - w) * Xk
      }
   }
   loss_bl = norm(x);
   loss_bl = deltanorm(xkp, x) / loss_bl;

}

// метод блочной релаксации: получение решения
template <typename T>
void MatrixSystem<T>::block_relaxation(std::ofstream& fs, T w)
{
   if (n % block_size != 0)
      throw 1;
   T loss = 0;
   for (int k = 0; k < n; k++)
      x[k] = 0;

   std::vector<T> buf(n);
   size_t j = 0;

   std::cout << w << std::endl;
   for (j = 0; j < max_iter; j++)
   {
      block_iteration(w, buf, loss);
      x.swap(buf);
      if (loss < eps)
      {
         std::cout << "End iter:" << j << std::endl;
         break;
      }
      if (j % 100 == 99)
      {
         std::cout << "iter: " << j + 1 << std::endl;
         std::cout << "loss: " << loss << std::endl;
      }
   }

   fs << w << "\t" << x[0] << "\t" << j + 1 << "\t" << loss << std::endl;
   for (j = 1; j < n; j++)
   {
      fs << "\t" << x[j] << std::endl;
   }

}

// подсчет числа обусловленностей для метода блочной релаксации
template <typename T>
T MatrixSystem<T>::num_bl_obusl()
{
   T obusl = 0;
   T r = 0;
   T sum = 0;
   
   multMV(x, buf);
   for (int i = 0; i < n; i++)
   {
      r = b[i] - buf[i];
      sum += r * r;
   }
   obusl = sqrt(sum);
   obusl = obusl / normb;
   sum = 0;
   for (int i = 0; i < n; i++)
   {
      r = x[i] - (i + 1);
      sum += r * r;
   }
   sum = sqrt(sum);
   obusl = (sum / sqrt(650.0)) / obusl;
   return obusl;
}

// норма разности двух векторов
template <typename T>
T MatrixSystem<T>::deltanorm(std::vector<T> x, std::vector<T> y)
{
   T norma = 0;
   for (size_t i = 0; i < n; i++)
      norma += (x[i] - y[i]) * (x[i] - y[i]);

   norma = sqrt(norma);
   return norma;

}

// умножение матрицы на вектор
template <typename T>
void MatrixSystem<T>::multMV(std::vector<T>& x, std::vector<T>& res)
{
   for (int i = 0; i < n; i++)
      res[i] = multVV(3, i, x);
}

// умножение блоков на вектор
template <typename T>
void MatrixSystem<T>::multBV(int blocknumber, std::vector<T>& xkp)
{
   int j = block_size * blocknumber;
   int y = j;
   for (int l = 0; l < block_size; l++, y++)
   {
      xkp[y] = 0;
   }
   if (j > 0)
   {
      xkp[j] += al[0][j - 1] * xkp[j - 1];
   }
   for (int l = 0; l < block_size; l++, j++)
   {
      //умножение нижнего треугольника матрицы
      if (j > m + 1)
      {
         xkp[j] += al[1][j - m - 2] * xkp[j - m - 2];
         if (j > m + k + 2)
         {
            xkp[j] += al[2][j - m - k - 3] * xkp[j - m - k - 3];
         }
      }

      // умножение верхнего треугольника матрицы
      if (j < n - m - 2)
      {
         xkp[j] += au[1][j] * xkp[j + m + 2];
         if (j < n - m - k - 3)
         {
            xkp[j] += au[2][j] * xkp[j + m + k + 3];
         }
      }

   }
   j--;
   if (j < n - 1)
   {
      xkp[j] += au[0][j] * xkp[j + 1];
   }

}
