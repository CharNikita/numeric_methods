#pragma once

#include <iomanip>

template<typename T>
class MatrixSystem
{
   int* ia, * ja, max_iter, n, flag;
   T* al, * au, * di;
   T* al_new, * au_new, * di_new;
   T* b;
   T* r, * z, * p;
   T* temp1, * temp2;
   T* MV;
   T skal1, skal2;
   T eps;

public:
   T* x;

   ~MatrixSystem();

   int getSize() { return n; }

   MatrixSystem(std::string file_kuslau, std::string file_gi, 
      std::string file_gj, std::string file_di, std::string file_ggu, 
      std::string file_ggl, std::string file_b, int flag);
   MatrixSystem(std::string file_params, std::string file_ia, int flag);

   void read_vector_T(std::string name, T* arr, int n);
   void read_vector_int(std::string name, int* arr, int n);
   void mult(T* MV, T* vec);
   void multU(T* MV, T* vec);
   void mult_tr(T* MV, T* vec);
   T scalar_mult(T* vec1, T* vec2);
   T norm(T* vec);
   void mult_pr(T* aa, T* y, T* b);
   void mult_obr(T* aa, T* y, T* b);
   void lu_sq();
   void diag();
   void msg_basic();
   void los_basic();
   void msg_sq();
   void los_sq();
   void msg_di();
   void los_di();
   void gilbert_generate(int n);
};