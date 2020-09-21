#include "MatrixSystem.cpp"

int main()
{
	std::string path = "input.txt";
	auto ms = MatrixSystem<double>(path);
	ms.ldu();
	ms.forward_pass();
	ms.central_pass();
	ms.backward_pass();
	ms.output();

   return 0;
}