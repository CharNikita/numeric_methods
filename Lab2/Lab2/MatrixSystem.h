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
	std::vector<T> buf;
	size_t n;
	size_t m;
	size_t k;
	size_t max_iter;
	double eps;
	size_t block_size;
	T normb;

	void readFromFile(std::string& path);
	MatrixSystem();

public:
	size_t get_size() { return n; }

	MatrixSystem(std::string& path);

	std::vector<T> get_x();
	void set_x(std::vector<T> x);
	T norm(std::vector<T>& x);
	T multVV(int flag, int i, std::vector<T>& x0);
	T jacobi(T w, std::vector<T>& x1, T& loss);
	T gauss_zeidel(T w, T& loss);
	T num_obusl(std::vector<T> x, T loss, T normxstar);
	void iteration(std::string& path, bool flag);
	void lu_factorization();
	void lu_solution(std::vector<T>& x0, int block_amt, T w);
	T num_bl_obusl();
	T deltanorm(std::vector<T> x, std::vector<T> y);
	void multMV(std::vector<T>& x, std::vector<T>& res);
	void multBV(int blocknumber, std::vector<T>& xkp);
	void block_iteration(T w, std::vector<T>& xkp, T& loss_bl);
	void block_relaxation(std::ofstream& path, T w);
};