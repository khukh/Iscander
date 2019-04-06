#include "pch.h"
#include "RGParam.h"
//#include "matrix.cpp"

RGParam::RGParam(double pitch, double yaw, double roll) {
	RGPar.resize(4);
	RGAngle.resize(3);
	RGAngle[0] = pitch;
	RGAngle[1] = yaw;
	RGAngle[2] = roll;

	RGPar[0] = cos(0.5*RGAngle[1])*cos(0.5*RGAngle[0])*cos(0.5*RGAngle[2]) - sin(0.5*RGAngle[1])*sin(0.5*RGAngle[0])*sin(0.5*RGAngle[2]);
	RGPar[1] = sin(0.5*RGAngle[1])*sin(0.5*RGAngle[0])*cos(0.5*RGAngle[2]) + cos(0.5*RGAngle[1])*cos(0.5*RGAngle[0])*sin(0.5*RGAngle[2]);
	RGPar[2] = sin(0.5*RGAngle[1])*cos(0.5*RGAngle[0])*cos(0.5*RGAngle[2]) + cos(0.5*RGAngle[1])*sin(0.5*RGAngle[0])*sin(0.5*RGAngle[2]);
	RGPar[3] = cos(0.5*RGAngle[1])*sin(0.5*RGAngle[0])*cos(0.5*RGAngle[2]) - sin(0.5*RGAngle[1])*cos(0.5*RGAngle[0])*sin(0.5*RGAngle[2]);

	
}


std::vector<double> RGParam::getRGPar() {
	std::vector <double> rgPar = RGPar;
	//?
	return RGPar;
}

RGParam::~RGParam() {}
