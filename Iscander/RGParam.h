#pragma once
#include "pch.h"
class RGParam {
public:
	RGParam(double pitch, double yaw, double roll);
	std::vector <double> RGPar;
	//std::vector <double> RGAngle;

	
	std::vector <double> getRGPar();
	~RGParam();
};

