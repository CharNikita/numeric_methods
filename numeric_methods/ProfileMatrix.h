#pragma once
#include <vector>
#include <fstream>
#include <string>


template <typename T>
class MatrixSystem {
	std::vector<T> di;
	std::vector<T> ia;
	std::vector<T> al;
	std::vector<T> au;
	size_t size;

    void readFromFile(std::string path);
	MatrixSystem();
	
public:
	size_t getSize() { return size };
	MatrixSystem(String path);
	~MatrixSystem();
};
