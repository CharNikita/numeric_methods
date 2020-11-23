#include "MatrixSystem.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>

using namespace std;

template<typename T>
void MatrixSystem<T>::read_vector_T(std::string name, T* arr, int n)
{
   name += ".txt";
   ifstream fin(name);
   for (int i = 0; i < n; i++)
      fin >> arr[i];

   fin.close();
}

template<typename T>
void MatrixSystem<T>::read_vector_int(std::string name, int* arr, int n)
{
   name += ".txt";
   ifstream fin(name);
   for (int i = 0; i < n; i++)
      fin >> arr[i];

   fin.close();
}


template <typename T>
MatrixSystem<T>::MatrixSystem(std::string file_kuslau, std::string file_gi, 
   std::string file_gj, std::string file_di, std::string file_ggu, 
   std::string file_ggl, std::string file_b, int flag)
{
   ifstream fin;
   fin.open(file_kuslau + ".txt");
   fin >> n >> max_iter >> eps;
   fin.close();

   this->flag = flag;
   b = new T[n];
   di = new T[n];
   x = new T[n];
   read_vector_T(file_di, di, n);
   read_vector_T(file_b, b, n);

   for (int i = 0; i < n; i++)
      x[i] = 0;

   ia = new int[n + 1];
   read_vector_int(file_gi, ia, n + 1);

   if (ia[0] == 1)
      for (int i = 0; i < n + 1; i++)
         ia[i]--;

   int m = ia[n];
   ja = new int[m];
   al = new T[m];
   au = new T[m];
   read_vector_int(file_gj, ja, m);

   if (ja[0] == 1)
      for (int i = 0; i < m; i++)
         ja[i]--;

   read_vector_T(file_ggl, al, m);
   read_vector_T(file_ggu, au, m);
   switch (flag)
   {
      case 1: //МСГ без предобуславливания
      {
         r = new T[n];
         z = new T[n];
         temp1 = new T[n];
         temp2 = new T[n];
         break;
      }
      case 2: //ЛОС без предобуславливания
      {
         r = new T[n];
         z = new T[n];
         p = new T[n];
         temp1 = new T[n];
         break;
      }
      case 3: //МСГ с лу-предобуславливанием
      {
         r = new T[n];
         z = new T[n];
         temp1 = new T[n];
         temp2 = new T[n];

         di_new = new T[n];
         int m = ia[n];
         al_new = new T[m];
         au_new = new T[m];

         break;
      }
      case 4: //ЛОС с лу-предобуславливанием
      {
         r = new T[n];
         z = new T[n];
         p = new T[n];
         temp1 = new T[n];
         temp2 = new T[n];

         di_new = new T[n];
         int m = ia[n];
         al_new = new T[m];
         au_new = new T[m];

         break;
      }
      case 5: //МСГ с di-предобуславливанием
      {
         r = new T[n];
         z = new T[n];
         temp1 = new T[n];
         temp2 = new T[n];

         di_new = new T[n];
         break;
      }
      case 6: //ЛОС с di-предобуславливанием
      {
         r = new T[n];
         z = new T[n];
         p = new T[n];
         temp1 = new T[n];
         temp2 = new T[n];

         di_new = new T[n];
         break;
      }
      default: break;
   }
}


template <typename T>
MatrixSystem<T>::MatrixSystem(std::string file_params, std::string file_ia, int flag)
{
   ifstream fin(file_params + ".txt");
   fin >> n >> eps >> max_iter;
   fin.close();

   this->flag = flag;
   b = new T[n];
   di = new T[n];
   x = new T[n];
   for (int i = 0; i < n; i++)
      x[i] = 0;
   ia = new int[n + 1];
   read_vector_int(file_ia, ia, n + 1);
   int m = ia[n];

   ja = new int[m];
   al = new T[m];
   au = new T[m];

   switch (flag)
   {
      case 1: //МСГ без предобуславливания
      {
         r = new T[n];
         z = new T[n];
         temp1 = new T[n];
         temp2 = new T[n];
         break;
      }
      case 2: //ЛОС без предобуславливания
      {
         r = new T[n];
         z = new T[n];
         p = new T[n];
         temp1 = new T[n];
         break;
      }
      case 3://МСГ с лу-предобуславливанием
      {
         r = new T[n];
         z = new T[n];
         temp1 = new T[n];
         temp2 = new T[n];

         di_new = new T[n];
         int m = ia[n];
         al_new = new T[m];
         au_new = new T[m];
         break;
      }
      case 4: //ЛОС с лу-предобуславливанием
      {
         r = new T[n];
         z = new T[n];
         p = new T[n];
         temp1 = new T[n];
         temp2 = new T[n];

         di_new = new T[n];
         int m = ia[n];
         al_new = new T[m];
         au_new = new T[m];
         break;
      }
      case 5: //МСГ с di-предобуславливанием
      {
         r = new T[n];
         z = new T[n];
         temp1 = new T[n];
         temp2 = new T[n];

         di_new = new T[n];
         break;
      }
      case 6: //ЛОС с di-предобуславливанием
      {
         r = new T[n];
         z = new T[n];
         p = new T[n];
         temp1 = new T[n];
         temp2 = new T[n];

         di_new = new T[n];
         break;
      }
      default:
         break;
   }
   gilbert_generate(n);
}

template <typename T>
MatrixSystem<T>::~MatrixSystem()
{
   delete x;
   delete b;
   delete al;
   delete au;
   delete ia;
   delete ja;
   delete di;

   switch (flag)
   {
      case 1: //МСГ без предобуславливания
      {
         delete r;
         delete z;
         delete temp1;
         delete temp2;
         break;
      }
      case 2: //ЛОС без предобуславливания
      {
         delete r;
         delete z;
         delete p;
         delete temp1;
         break;
      }
      case 3: //МСГ с лу-предобуславливанием
      {
         delete r;
         delete z;
         delete temp1;
         delete temp2;

         delete di_new;
         delete al_new;
         delete au_new;
         break;
      }
      case 4: //ЛОС с лу-предобуславливанием
      {
         delete r;
         delete z;
         delete p;
         delete temp1;
         delete temp2;

         delete di_new;
         delete al_new;
         delete au_new;
         break;
      }
      case 5: //МСГ с di-предобуславливанием
      {
         delete r;
         delete z;
         delete temp1;
         delete temp2;
         delete di_new;
         break;
      }
      case 6: //ЛОС с di-предобуславливанием
      {
         delete r;
         delete z;
         delete p;
         delete temp1;
         delete temp2;
         delete di_new;
         break;
      }
      default:
         break;
   }
}

//умножение матрицы в разреженном на вектор
template <typename T>
void MatrixSystem<T>::mult(T* MV, T* vec)
{
   for (int i = 0; i < n; i++) 
   {
      int k0 = ia[i];
      int k1 = ia[i + 1];
      MV[i] = di[i] * vec[i];

      for (int k = k0; k < k1; k++) 
      {
         int j = ja[k];
         MV[i] += vec[j] * al[k];
         MV[j] += vec[i] * au[k];
      }
   }
}

template <typename T>
void MatrixSystem<T>::multU(T* MV, T* vec)
{
   for (int i = 0; i < n; i++) 
   {
      int k0 = ia[i];
      int k1 = ia[i + 1];
      MV[i] = di[i] * vec[i];

      for (int k = k0; k < k1; k++) 
      {
         int j = ja[k];
         //MV[i] += vec[j] * al[k];
         MV[j] += vec[i] * au_new[k];
      }
   }
}

template <typename T>
void MatrixSystem<T>::mult_tr(T* MV, T* vec)
{
   for (int i = 0; i < n; i++) 
   {
      int k0 = ia[i];
      int k1 = ia[i + 1];
      MV[i] = di[i] * vec[i];

      for (int k = k0; k < k1; k++) 
      {
         int j = ja[k];
         MV[i] += vec[j] * au[k];
         MV[j] += vec[i] * al[k];
      }
   }
}

template <typename T>
T MatrixSystem<T>::scalar_mult(T* vec1, T* vec2)
{
   T s = 0;
   for (int i = 0; i < n; i++) 
      s += vec1[i] * vec2[i];

   return s;
}

template <typename T>
T MatrixSystem<T>::norm(T* vec)
{
   T sum = 0;
   for (int i = 0; i < n; i++)
      sum += vec[i] * vec[i];

   return sqrt(sum);
}

template <typename T>
void MatrixSystem<T>::mult_pr(T* aa, T* y, T* b)
{

   for (int i = 0; i < n; i++) 
   {
      T s = 0; //переменные суммирования

      int i0 = ia[i];//индекс 1го элемента в iтой строке
      int i1 = ia[i + 1];

      for (int k = i0; k < i1; k++) 
      {
         int j = ja[k];
         s += y[j] * aa[k];
      }

      y[i] = (b[i] - s) / di_new[i];
   }
}

template <typename T>
void MatrixSystem<T>::mult_obr(T* aa, T* y, T* b)
{
   for (int i = 0; i < n; i++)
      y[i] = b[i];

   for (int i = n - 1; i >= 0; i--) 
   {
      int i0 = ia[i];//индекс 1го элемента в iтой строке
      int i1 = ia[i + 1];

      y[i] /= di_new[i];

      for (int k = i1 - 1; k >= i0; k--) 
      {
         int j = ja[k];
         y[j] -= y[i] * aa[k];
      }
   }
}

template <typename T>
void MatrixSystem<T>::lu_sq()
{
   //копирование-инициализация
   for (int i = 0; i < n; i++)
      di_new[i] = di[i];

   for (int i = 0; i < ia[n]; i++) 
   {
      al_new[i] = al[i];
      au_new[i] = au[i];
   }

   for (int i = 0; i < n; i++) 
   {
      T sum = 0; //переменные суммирования

      int i0 = ia[i];
      int i1 = ia[i + 1];

      for (int k = i0; k < i1; k++) 
      {
         int j = ja[k];
         T sl = 0, su = 0;
         int j0 = ia[j];
         int j1 = ia[j + 1];
         int ki = i0;
         int kj = j0;

         for (; ki < k && kj < j1;) 
         {
            int jl = ja[ki];
            int ju = ja[kj];

            if (jl == ju) 
            {
               sl += au_new[kj] * al_new[ki];
               su += al_new[kj] * au_new[ki];
               ki++; kj++;
            }

            else if (jl < ju) ki++;
            else kj++;
         }

         au_new[k] = (au_new[k] - su) / di_new[j];
         al_new[k] = (al_new[k] - sl) / di_new[j];
         sum += au_new[k] * al_new[k];
      }

      di_new[i] = sqrt(di_new[i] - sum);
   }
}

template <typename T>
void MatrixSystem<T>::diag()
{
   for (int i = 0; i < n; i++)
      di_new[i] = (double)1 / sqrt(di[i]);
}

template <typename T>
void MatrixSystem<T>::msg_basic()
{
   //	инициализация
   mult(temp1, x);
   for (int i = 0; i < n; i++)
      temp1[i] = b[i] - temp1[i];

   mult_tr(r, temp1);

   for (int i = 0; i < n; i++)
      z[i] = r[i];

   skal1 = scalar_mult(r, r);
   T norm_f = norm(b);

   //iteration
   T nev = norm(r) / norm_f;
   for (int k = 0; k < max_iter && nev > eps; k++) 
   {
      cout << k + 1 << ": ";
      mult(temp2, z);

      mult_tr(temp1, temp2);
      skal2 = scalar_mult(temp1, z);

      T alpha = skal1 / skal2;
      for (int i = 0; i < n; i++) 
      {
         x[i] += alpha * z[i];
         r[i] -= alpha * temp1[i];
      }

      skal2 = scalar_mult(r, r);

      T beta = skal2 / skal1;
      for (int i = 0; i < n; i++)
         z[i] = r[i] + beta * z[i];

      skal1 = skal2;

      nev = norm(r) / norm_f;
      //nev = skal_mult(r, r) / norm_f;
//	   for (int i = 0; i < n; i++)
//		   cout << x[i] << " ";
      cout << "nev=" << nev << endl;
   }
}

template <typename T>
void MatrixSystem<T>::los_basic() {
   //	инициализация
   mult(temp1, x);
   for (int i = 0; i < n; i++) 
   {
      r[i] = b[i] - temp1[i];
      z[i] = r[i];
   }
   mult(p, z);

   //iteration
   T nev = scalar_mult(r, r);
   for (int k = 0; k < max_iter && abs(nev) > eps; k++) 
   {
      cout << k + 1 << ": ";
      skal1 = scalar_mult(p, r);
      skal2 = scalar_mult(p, p);

      T alpha = skal1 / skal2;

      for (int i = 0; i < n; i++) 
      {
         x[i] += alpha * z[i];
         r[i] -= alpha * p[i];
      }

      mult(temp1, r);

      skal1 = scalar_mult(p, temp1);

      T beta = -skal1 / skal2;
      for (int i = 0; i < n; i++) 
      {
         z[i] = r[i] + beta * z[i];
         p[i] = temp1[i] + beta * p[i];
      }
      //mult(p, z);

      //nev = nev - alpha * alpha*skal2;

      nev = scalar_mult(r, r);

      cout << "nev=" << nev << endl;
   }
};

template <typename T>
void MatrixSystem<T>::msg_sq()
{
   lu_sq();

   //	инициализация
   mult(temp1, x);
   for (int i = 0; i < n; i++)
      temp1[i] = b[i] - temp1[i];

   mult_pr(al_new, temp2, temp1);
   mult_obr(al_new, temp1, temp2);
   mult_tr(temp2, temp1);
   mult_pr(au_new, r, temp2);

   for (int i = 0; i < n; i++)
      z[i] = r[i];

   multU(temp1, x);
   for (int i = 0; i < n; i++)
      x[i] = temp1[i];


   skal1 = scalar_mult(r, r);
   T norm_f = norm(b);

   //iteration
   T nev = norm(r) / norm_f;

   for (int k = 0; k < max_iter && nev > eps; k++) 
   {
      cout << k + 1 << ": ";
      mult_obr(au_new, temp2, z);
      mult(temp1, temp2);
      mult_pr(al_new, temp2, temp1);
      mult_obr(al_new, temp1, temp2);
      mult_tr(temp2, temp1);
      mult_pr(au_new, temp1, temp2);

      skal2 = scalar_mult(temp1, z);

      T alpha = skal1 / skal2;
      for (int i = 0; i < n; i++) 
      {
         x[i] += alpha * z[i];
         r[i] -= alpha * temp1[i];
      }

      skal2 = scalar_mult(r, r);

      T beta = skal2 / skal1;
      for (int i = 0; i < n; i++)
         z[i] = r[i] + beta * z[i];

      skal1 = skal2;

      nev = norm(r) / norm_f;
      cout << "nev=" << nev << endl;
   }

   for (int i = 0; i < n; i++)
      temp1[i] = x[i];

   mult_obr(au_new, x, temp1);
}

template <typename T>
void MatrixSystem<T>::msg_di() {

   diag();

   //	инициализация
   mult(temp1, x);
   for (int i = 0; i < n; i++) 
   {
      temp1[i] = b[i] - temp1[i];
      temp1[i] /= di[i];
   }

   //mult_di(temp2,temp1);
  // mult_di(temp1, temp2);
   mult_tr(temp2, temp1);
   //mult_di(r, temp2);

   for (int i = 0; i < n; i++) 
   {
      r[i] = di_new[i] * temp2[i];
      z[i] = r[i];
      x[i] = x[i] * (double)sqrt(di[i]);
   }

   skal1 = scalar_mult(r, r);
   T norm_f = norm(b);

   //iteration
   T nev = norm(r) / norm_f;

   for (int k = 0; k < max_iter && nev > eps; k++) 
   {
      cout << k + 1 << ": ";

      for (int i = 0; i < n; i++)
         temp2[i] = z[i] * di_new[i];

      mult(temp1, temp2);

      for (int i = 0; i < n; i++)
         temp1[i] = temp1[i] / di[i];

      mult_tr(temp2, temp1);
      for (int i = 0; i < n; i++)
         temp1[i] = temp2[i] * di_new[i];

      skal2 = scalar_mult(temp1, z);

      T alpha = skal1 / skal2;
      for (int i = 0; i < n; i++) 
      {
         x[i] += alpha * z[i];
         r[i] -= alpha * temp1[i];
      }

      skal2 = scalar_mult(r, r);

      T beta = skal2 / skal1;
      for (int i = 0; i < n; i++) 
         z[i] = r[i] + beta * z[i];

      skal1 = skal2;

      nev = norm(r) / norm_f;

      cout << "nev=" << nev << endl;
   }

   for (int i = 0; i < n; i++)
      x[i] = x[i] * di_new[i];

}

template <typename T>
void MatrixSystem<T>::los_sq() {

   lu_sq();

   //	инициализация
   mult(temp1, x);
   for (int i = 0; i < n; i++)
      temp2[i] = b[i] - temp1[i];

   mult_pr(al_new, r, temp2);

   mult_obr(au_new, z, r);

   mult(temp1, z);
   mult_pr(al_new, p, temp1);

   //iteration
   T nev = scalar_mult(r, r);
   for (int k = 0; k < max_iter && nev > eps; k++) 
   {
      cout << k + 1 << ": ";
      skal1 = scalar_mult(p, r);
      skal2 = scalar_mult(p, p);

      T alpha = skal1 / skal2;
      for (int i = 0; i < n; i++) 
      {
         x[i] += alpha * z[i];
         r[i] -= alpha * p[i];
      }
      mult_obr(au_new, temp1, r);
      mult(temp2, temp1);
      mult_pr(al_new, temp1, temp2);

      skal1 = scalar_mult(p, temp1);

      T beta = -skal1 / skal2;

      mult_obr(au_new, temp2, r);

      for (int i = 0; i < n; i++)
         z[i] = temp2[i] + beta * z[i];

      for (int i = 0; i < n; i++)
         p[i] = temp1[i] + beta * p[i];

      //   nev = nev - alpha * alpha*skal2;
      nev = scalar_mult(r, r);
      cout << "nev=" << nev << endl;
   }
};

template <typename T>
void MatrixSystem<T>::los_di()
{
   diag();

   //	инициализация
   mult(temp1, x);
   for (int i = 0; i < n; i++) 
   {
      temp1[i] = b[i] - temp1[i];
      r[i] = temp1[i] * di_new[i];
      z[i] = di_new[i] * r[i];
   }

   mult(temp1, z);

   for (int i = 0; i < n; i++)
      p[i] = di_new[i] * temp1[i];

   //iteration
   T nev = scalar_mult(r, r);
   for (int k = 0; k < max_iter && nev > eps; k++) 
   {
      cout << k + 1 << ": ";
      skal1 = scalar_mult(p, r);
      skal2 = scalar_mult(p, p);

      T alpha = skal1 / skal2;

      for (int i = 0; i < n; i++) 
      {
         x[i] += alpha * z[i];
         r[i] -= alpha * p[i];
         temp1[i] = di_new[i] * r[i];
      }

      mult(temp2, temp1);

      for (int i = 0; i < n; i++)
         temp1[i] = di_new[i] * temp2[i];

      skal1 = scalar_mult(p, temp1);

      T beta = -skal1 / skal2;

      for (int i = 0; i < n; i++) 
      {
         z[i] = di_new[i] * r[i] + beta * z[i];
         p[i] = temp1[i] + beta * p[i];
      }

      nev = scalar_mult(r, r);
      //nev = nev - alpha * alpha*skal2;
      cout << "nev=" << nev << endl;
   }
};

template <typename T>
void MatrixSystem<T>::gilbert_generate(int n)
{

   for (int i = 0; i < n; i++)
      b[i] = 0;

   for (int i = 0; i < n; i++) 
   {
      int k = ia[i];
      for (int j = 0; j < i; j++, k++) 
      {
         ja[k] = j;
         al[k] = 1.0 / (i + j + 1);
         au[k] = 1.0 / (i + j + 1);
         //b[i] += au[k] * x_teor[j];
         b[i] += au[k] * (j + 1);
         b[j] += al[k] * (i + 1);
      }

      di[i] = 1.0 / (i + i + 1);
      b[i] += di[i] * (i + 1);
   }

   ofstream gout;
   gout.open("di.txt");
   gout << setprecision(17);
   for (int i = 0; i < n; i++)
      gout << di[i] << " ";
   gout << endl;
   gout.close();

   gout.open("gj.txt");
   gout << setprecision(17);
   for (int i = 1; i < ia[n]; i++)
      gout << ja[i] + 1 << " ";
   gout << endl;
   gout.close();

   gout.open("ggl.txt");
   gout << setprecision(17);
   for (int i = 1; i < ia[n]; i++)
      gout << al[i] << " ";
   gout << endl;
   gout.close();

   gout.open("ggu.txt");
   gout << setprecision(17);
   for (int i = 1; i < ia[n]; i++)
      gout << au[i] << " ";
   gout << endl;
   gout.close();

   gout.open("pr.txt");
   gout << setprecision(17);
   for (int i = 0; i < n; i++)
      gout << b[i] << " ";
   gout << endl;
   gout.close();
}