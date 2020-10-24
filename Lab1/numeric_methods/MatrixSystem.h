#pragma once

#include <vector>
#include <string>

template <typename T>
class MatrixSystem
{
	std::vector<T> di;
	std::vector<T> ia;
	std::vector<T> al;
	std::vector<T> au;
	std::vector<T> b;
	size_t size;

   void readFromFile(std::string path);
	MatrixSystem();
	
public:
	size_t getSize() { return size; }

	MatrixSystem(std::string path);

	~MatrixSystem();

	void ldu();
	void forward_pass();
	void central_pass();
	void backward_pass();
	void output();
};
