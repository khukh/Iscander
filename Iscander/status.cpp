
#include "pch.h"

#define _USE_MATH_DEFINES
#include "status.h"

#include "cmath"
#include <fstream>


status::status(): Rot(PITCH0, YAW0, ROLL0) {
	parametr.resize(14);
	ForcePr.resize(3);
	ForcePrG.resize(3);
	ForceG.resize(3);
	Torque.resize(3);
	//ADkoef.resize(6);
	//скорости
	parametr[0] = V0 * cos(PITCH0) * cos(YAW0);
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
	//параметры Родрига-Гамильтона
	std::vector <double> a = Rot.RG.getRGPar();
	parametr[9] = a[0];
	parametr[10] = a[1];
	parametr[11] = a[2];
	parametr[12] = a[3];

	parametr[13] = 0;
}

status::~status() {}

//значения производных
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
	prir[11] = 0.5*(parametr[6] * parametr[12] + parametr[7] * parametr[9] - parametr[8] * parametr[10]);
	prir[12] = 0.5*(-parametr[6] * parametr[11] + parametr[7] * parametr[10] + parametr[8] * parametr[9]);
	//dt/dt
	prir[13] = 1;

	return(prir);

}

void status::nonIntegr() {
	Rot.RG.setRGPar(parametr[9], parametr[10], parametr[11], parametr[12]);
	Rot.RG.norm();
	Rot.fromRGtoAngles();

	std::vector <double> v(3);
	std::vector <double> vg {parametr[0],parametr[1],parametr[2]};
	Rot.fromRGtoMatrixT();
	v = Rot.A*vg;

	double vFullsq = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	density = GOST4401.roFunc(parametr[4]);
	double ah = atan2(v[1], v[0]);
	alpha = - ah;
	betta = asin(v[2] / sqrt(vFullsq));
	double alphaSpace = sqrt(alpha*alpha + betta * betta);
	double mach = sqrt(vFullsq) / GOST4401.aFunc(parametr[4]);

	double deltaT = K_T1 * (parametr[8]*cos(Rot.Angles[2]) + parametr[7]*sin(Rot.Angles[2]));
	double deltaN = K_N1 * (parametr[7] * cos(Rot.Angles[2]) - parametr[8] * sin(Rot.Angles[2]))/cos(Rot.Angles[0]);
	double deltaE = K_E1 * Rot.Angles[2] + K_E2 * (parametr[6] - tan(Rot.Angles[0])*(parametr[7]*cos(Rot.Angles[2])-parametr[8]*sin(Rot.Angles[2])));

	double q = density * vFullsq / 2;
	ForcePr[0] = -Cx(mach, alphaSpace) * density * vFullsq * S_M / 2;
	ForcePr[1] = (CyAl(mach, alpha)*alpha +CyDeltaT(mach, alpha)*deltaT) * density * vFullsq * S_M / 2;
	ForcePr[2] = (CzBetta(mach, betta)*betta + CzDeltaN(mach, betta)*deltaN) * density * vFullsq * S_M / 2;
	Rot.fromRGtoMatrix();
	ForcePrG = Rot.A * ForcePr;
	g = 9.80665/*PI0 / ((R_EARTH_G + parametr[4])*(R_EARTH_G + parametr[4]))*/;
	ForceG[0] = 0;
	ForceG[1] = -M * g;
	ForceG[2] = 0;

	//TODO: вставить моменты!
	//double q = L * density * vFullsq * S_M / 2;
	double kADdevelopment = 0.0;
	double kMVr = 0;
	double P = 0; //тяга?
	double Mstab = 0.5 * D_M * P * deltaE;
	Torque[0] = (MxOmegaX(mach, alphaSpace) * parametr[6] * L / sqrt(vFullsq) + MxDelta(mach, alpha, q*S_M*L) * deltaE)* density * vFullsq * S_M * L / 2 + Mstab;
	Torque[1] = (MyOmegaY(mach, betta) * parametr[7] * L / sqrt(vFullsq) + MyBetta(mach, betta) * betta + MyDelta(mach, betta, q*S_M*L) * deltaN) * density * vFullsq * S_M * L / 2;
	Torque[2] = (MzOmegaZ(mach, alpha) * parametr[8] * L / sqrt(vFullsq) + MzAlpha(mach, alpha) * alpha + MzDelta(mach, alpha, q*S_M*L) *deltaT) * density * vFullsq * S_M * L / 2;

	


}

//вывод параметров
void status::printParam(std::ofstream &fout) {
	fout << parametr[13] << '\t' << std::scientific;
	for (int i = 0; i <13; i++) {
		fout << parametr[i] << '\t';
		
	}
	
	fout << alpha * 180 / PI << '\t' << betta << '\t';

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
		fout << Rot.Angles[i]*180/PI << '\t';
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

	Rot = right.Rot;

	alpha = right.alpha;
	betta = right.betta;

	g = right.g;
	//ADkoef = right.ADkoef;
	
	density = right.density;

	return *this;
}


