#pragma once

#include <vector>
#include <string>

template <typename T>
class MatrixSystem
{
	std::vector<T> di;
	std::vector<std::vector<T>> al;
	std::vector<std::vector<T>> au;
	std::vector<T> b;
	std::vector<T> x;
	size_t n;
	size_t m;
	size_t k;
	size_t max_iter;
	double eps;
	size_t block_size;
	T normb;

	void readFromFile(std::string path);
	MatrixSystem();

public:
	size_t getSize() { return n; }

	MatrixSystem(std::string path);

	T norm(std::vector<T> x);
	T multVV(int flag, int i, std::vector<T> x0);
	T jacobi(T w, std::vector<T>& x1, T& loss);
	T gaussZeidel(T w, std::vector<T>& x1, T& loss);
	T num_obusl(std::vector<T> x, T loss, T normxstar);
	void iteration(std::string path, bool flag);
	void lu();
	void output();
};