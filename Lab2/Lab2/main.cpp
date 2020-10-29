#include "MatrixSystem.cpp"
//#include "GaussMatrixSystem.cpp"

#include <iostream>

int main()
{
	std::string path = "input.txt";
	auto ms = MatrixSystem<double>(path);

	ms.iteration("output.txt", true);

   //double w = 0;
   //double obusl;
   //for (int i = 120; i < 140; ++i)
   //{
   //   w = i / 100.0;
   //   ms.iteration("output.txt", false);
   //}

	//try
	//{
	//	ms.ldu();
	//}
	//catch (std::runtime_error e)
	//{
	//	std::cout << e.what();
	//	exit(1);
	//}

	//ms.output();

	return 0;
}