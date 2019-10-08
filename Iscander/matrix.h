#pragma once
#include "pch.h"
class matrix {
public:
	matrix(double pitch, double yaw, double roll);

	std::vector<std::vector<double>> matr;
	friend std::vector <double> operator*(matrix A, std::vector <double> b);
	//std::vector <double> operator * ( matrix &A, std::vector <double> &b);
	~matrix();
};

