
#include "pch.h"

#define _USE_MATH_DEFINES
#include "status.h"

#include "cmath"
#include <fstream>


status::status() {
	parametr.resize(13);
	
	//скорости
	parametr[0] = V0 * cos(PITCH0) * cos(YAW0);;
	parametr[1] = V0 * sin(PITCH0);
	parametr[2] = V0 * cos(PITCH0) * sin(YAW0);
	//кординаты
	parametr[3] = X0;
	parametr[4] = Y0;
	parametr[5] = Z0;
	//угловые скорости
	parametr[6] = 0;
	parametr[7] = 0;
	parametr[8] = 0;
	//параметры –одрига-√амильтона
	std::vector <double> a = rg.getRGPar();
	parametr[9] = 0;
	parametr[10] = 0;
	parametr[11] = 0;
	parametr[12] = 0;
}


status::~status() {}


//значени€ производных
std::vector <double> status::rightPart(int n) {
	std::vector <double> prir(10);
	
	mu = MU_SUN / pow(r, 3);
	prir[0] = -mu * parametr[3] + u[0];
	prir[1] = -mu * parametr[4] + u[1];
	prir[2] = -mu * parametr[5] + u[2];

	prir[3] = parametr[0];
	prir[4] = parametr[1];
	prir[5] = parametr[2];

	prir[6] = 1;
	prir[7] = u[0] * u[0];
	prir[8] = u[1] * u[1];
	prir[9] = u[2] * u[2];

	return(prir);

}



void status::nonIntegr(int n) {
	r = pow((pow(parametr[3], 2) + pow((parametr[4]), 2) + pow(parametr[5], 2)), 0.5);
	setUcontrol(n);
}


//вывод параметров
void status::printParam(std::ofstream &fout) {
	fout << '\t' << std::scientific << parametr[6] / 3600 / 24 << '\t' << parametr[0] << '\t' << parametr[1] << '\t' << parametr[2] << '\t' << parametr[3] << '\t' << parametr[4] << '\t' << parametr[5] << '\t' << std::scientific << u[0] << '\t' << std::scientific << u[1] << '\t' << std::scientific << u[2] << '\t' << std::scientific << parametr[7] << '\t' << std::scientific << parametr[8] << '\t' << std::scientific << parametr[9] << '\t' << r << '\t' << u[0] - parametr[0] << '\t' << u[1] - parametr[1] << '\t' << u[2] - parametr[2] << '\t' << mu << '\t' << MU_SUN / pow(r, 3) << '\n';

}
void status::setParam(std::vector <double> b) {
	parametr = b;
}

void status::setPSI(std::vector<std::vector<double>> psixNew, std::vector<std::vector<double>> psivNew) {
	psiX = psixNew;
	psiV = psivNew;

}

void status::setUcontrol(int n) {
	for (int i = 0; i < 3; i++) {
		u[i] = Ucontrol(n, parametr[6] - n * (SUM_T / N), psiX[n][i], psiV[n][i]);
	}

}
std::vector <double> status::getParam() {
	std::vector <double> get = parametr;
	return(get);
}
status& status::operator=(const status& right) {
	//проверка на самоприсваивание
	if (this == &right) {
		return *this;
	}
	parametr = right.parametr;
	r = right.r;
	u = right.u;
	psiX = right.psiX;
	psiV = right.psiV;

	return *this;
}


