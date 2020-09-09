#include "MatrixSystem.cpp"

int main()
{
	std::string path = "input.txt";
	auto ms = MatrixSystem<int>(path);
	ms.ldu();
   return 0;
}