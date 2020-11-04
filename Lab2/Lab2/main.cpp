#include "MatrixSystem.cpp"

#include <iostream>

int main()
{
	std::string in = "input.txt";
	std::string out = "output.txt";
	auto ms = MatrixSystem<double>(in);

	//ms.iteration(out, 2);

   //double w = 0;
   double obusl;
	ms.lu_factorization();
   for (int i = 160; i < 170; ++i)
   {
		double w = i / 100.0;
      ms.block_relaxation(out, w);
   }

	return 0;
}