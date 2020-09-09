#include "MatrixSystem.cpp"

int main()
{
	std::string path = "input.txt";
	auto ms = new MatrixSystem<int>(path);
	ms->ldu();
   return 0;
}