#include "pch.h"
#include "matrix.h"


matrix::matrix(double pitch, double yaw, double roll) { //из связанной в стартовую
	matr.assign(3, std::vector<double>(3));

	matr[0][0] = cos(pitch)*cos(yaw);
	matr[0][1] = -sin(pitch)*cos(yaw)*cos(roll) + sin(yaw)*sin(roll);
	matr[0][2] = sin(pitch)*cos(yaw)*sin(roll)+sin(yaw)*cos(roll);

	matr[0][0] = sin(pitch);
	matr[0][1] = cos(pitch)*cos(roll);
	matr[0][2] = -cos(pitch)*sin(roll);

	matr[0][0] = -cos(pitch)*sin(yaw);
	matr[0][1] = sin(pitch)*sin(yaw)*cos(roll)+cos(yaw)*sin(roll);
	matr[0][2] = -sin(pitch)*sin(yaw)*sin(roll) + cos(yaw)*cos(roll);
}




matrix::~matrix() {}

std::vector<double> operator*(matrix A, std::vector<double> b) {
	std::vector<double> res(3);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			res[i] += A.matr[i][j] * b[j];
		}

	}
	return res;
}
