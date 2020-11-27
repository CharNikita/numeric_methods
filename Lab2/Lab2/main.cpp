#include "MatrixSystem.cpp"

#include <iostream>

int main()
{
   std::string in = "input.txt";
   std::string out = "output.txt";
   auto ms = MatrixSystem<double>(in);

   //ms.iteration(out, true);

   std::ofstream fs;
   fs.open(out);
   fs.imbue(std::locale("Russian"));
   fs.precision(17);

   ms.lu_factorization();
   for (int i = 160; i < 170; ++i)
   {
      double w = i / 100.0;
      ms.block_relaxation(fs, w);

      auto ms_nev = MatrixSystem<double>(in);
      ms_nev.set_x(ms.get_x());
      double obusl = ms_nev.num_bl_obusl();
      fs << "Obusl: " << obusl << std::endl << std::endl;
   }

   fs.close();

   return 0;
}