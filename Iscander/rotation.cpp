#include "pch.h"
#include "rotation.h"


rotation::rotation(double pitch, double yaw, double roll) :A(pitch, yaw, roll), RG(pitch, yaw, roll) {
	Angles.resize(3);
	Angles[0] = pitch;
	Angles[0] = yaw;
	Angles[0] = roll;
}



void rotation::fromRGtoAngles() {
	Angles[0] = asin(2 * (RG.RGPar[0] * RG.RGPar[3] + RG.RGPar[1] * RG.RGPar[2]));
	Angles[1] = atan2(2 * (RG.RGPar[0] * RG.RGPar[1] - RG.RGPar[3] * RG.RGPar[2]), (pow(RG.RGPar[0], 2) + pow(RG.RGPar[2], 2) - pow(RG.RGPar[3], 2) - pow(RG.RGPar[1], 2)));
	Angles[2] = atan2(2 * (RG.RGPar[0] * RG.RGPar[2] - RG.RGPar[1] * RG.RGPar[3]), (pow(RG.RGPar[0], 2) + pow(RG.RGPar[1], 2) - pow(RG.RGPar[3], 2) - pow(RG.RGPar[2], 2)));
}

void rotation::fromRGtoMatrix() {
	A.matr[0][0] = (pow(RG.RGPar[0], 2) + pow(RG.RGPar[1], 2) - pow(RG.RGPar[2], 2) - pow(RG.RGPar[3], 2));
	A.matr[0][1] = 2*(RG.RGPar[0] * RG.RGPar[3] + RG.RGPar[1] * RG.RGPar[2]);
	A.matr[0][2] = 2 * (-RG.RGPar[0] * RG.RGPar[2] + RG.RGPar[1] * RG.RGPar[3]);

	A.matr[1][0] = 2 * (-RG.RGPar[0] * RG.RGPar[3] + RG.RGPar[1] * RG.RGPar[2]);
	A.matr[1][1] = (pow(RG.RGPar[0], 2) + pow(RG.RGPar[2], 2) - pow(RG.RGPar[3], 2) - pow(RG.RGPar[1], 2));
	A.matr[1][2] = 2 * (RG.RGPar[0] * RG.RGPar[1] + RG.RGPar[3] * RG.RGPar[2]);

	A.matr[2][0] = 2 * (RG.RGPar[0] * RG.RGPar[2] + RG.RGPar[1] * RG.RGPar[3]);
	A.matr[2][1] = 2 * (-RG.RGPar[0] * RG.RGPar[1] + RG.RGPar[2] * RG.RGPar[3]);
	A.matr[2][2] = (pow(RG.RGPar[0], 2) + pow(RG.RGPar[3], 2) - pow(RG.RGPar[1], 2) - pow(RG.RGPar[2], 2));
}


rotation::~rotation() {}
