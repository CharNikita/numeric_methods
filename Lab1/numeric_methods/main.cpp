#include "MatrixSystem.cpp"
#include "GaussMatrixSystem.cpp"

#include <iostream>

int main()
{
	std::string path = "input.txt";
	auto ms = MatrixSystem<double>(path);
	try
	{
		ms.ldu();
	}
	catch(std::runtime_error e)
	{
		std::cout << e.what();
		exit(1);
	}

	ms.forward_pass();
	ms.central_pass();
	ms.backward_pass();
	ms.output();

	/*auto ms = GaussMatrixSystem<double>(path);
	try
	{
		ms.gauss();
	}
	catch (std::runtime_error e)
	{
		std::cout << e.what();
		exit(1);
	}
	ms.output();*/

   return 0;
}