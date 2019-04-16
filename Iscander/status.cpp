
#include "pch.h"

#define _USE_MATH_DEFINES
#include "status.h"

#include "cmath"
#include <fstream>


status::status(): A(PITCH0, YAW0, ROLL0), Rg(PITCH0, YAW0, ROLL0) {
	parametr.resize(14);
	ForcePr.resize(3);
	ForcePrG.resize(3);
	ForceG.resize(3);
	Torque.resize(3);
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
	std::vector <double> a = Rg.getRGPar();
	parametr[9] = a[0];
	parametr[10] = a[1];
	parametr[11] = a[2];
	parametr[12] = a[3];

	parametr[13] = 0;
}

status::~status() {}

//значени€ производных
std::vector <double> status::rightPart() {
	std::vector <double> prir(14);
	//dV/dt
	prir[0] = (ForcePrG[0] + ForceG[0]) / M;
	prir[1] = (ForcePrG[1] + ForceG[1]) / M;
	prir[2] = (ForcePrG[2] + ForceG[2]) / M;
	//dX/dt
	prir[3] = parametr[0];
	prir[4] = parametr[1];
	prir[5] = parametr[2];
	//dW/dt
	prir[6] = Torque[0] / I_X - (I_Z - I_Y) / I_X * parametr[7] * parametr[8];
	prir[7] = Torque[1] / I_Y - (I_X - I_Z) / I_Y * parametr[6] * parametr[8];
	prir[8] = Torque[2] / I_Z - (I_Y - I_X) / I_Z * parametr[6] * parametr[7];
	//dPRG/dt
	prir[9] = -0.5*(parametr[6] * parametr[10] + parametr[7] * parametr[11] + parametr[8] * parametr[12]);
	prir[10] = 0.5*(parametr[6] * parametr[9] - parametr[7] * parametr[12] + parametr[8] * parametr[11]);
	prir[9] = 0.5*(parametr[6] * parametr[12] + parametr[7] * parametr[9] - parametr[8] * parametr[10]);
	prir[9] = 0.5*(-parametr[6] * parametr[11] + parametr[7] * parametr[10] + parametr[8] * parametr[9]);
	//dt/dt
	prir[13] = 1;

	return(prir);

}

void status::nonIntegr() {
	double vFullsq = parametr[0] * parametr[0] + parametr[1] * parametr[1] + parametr[2] * parametr[2];
	ForcePr[0] = -Cx * density * vFullsq * S_M / 2;
	ForcePr[1] = Cy * density * vFullsq * S_M / 2;
	ForcePr[2] = Cz * density * vFullsq * S_M / 2;

	ForcePrG = A * ForcePr;

	ForceG[0] = 0;
	ForceG[1] = -M * g;
	ForceG[2] = 0;

	alpha = -atan2(parametr[1], parametr[0]);
	betta = asin(parametr[2] / sqrt(vFullsq));
}

//вывод параметров
void status::printParam(std::ofstream &fout) {
	fout << parametr[13] << '\t' << std::scientific;
	for (int i = 0; i <13; i++) {
		fout << parametr[i] << '\t';
		
	}
	
	fout << alpha << '\t' << betta << '\t';

	for (int i = 0; i < 3; i++) {
		fout << ForceG[i] << '\t';
	}
	for (int i = 0; i < 3; i++) {
		fout << ForcePr[i] << '\t';
	}
	for (int i = 0; i < 3; i++) {
		fout << Torque[i] << '\t';
	}
	for (int i = 0; i < 3; i++) {
		fout << Rg.RGAngle[i] << '\t';
	}
	fout << '\n';

}
void status::setParam(std::vector <double> b) {
	parametr = b;
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

	ForcePr = right.ForcePr;
	ForcePrG = right.ForcePrG;  
	ForceG = right.ForceG;  
	Torque = right.Torque;  

	A = right.A;
	Rg = right.Rg;

	alpha = right.alpha;
	betta = right.betta;

	g = right.g;
	Cx = right.Cx;
	Cy = right.Cy;
	Cz = right.Cz;
	density = right.density;

	return *this;
}


