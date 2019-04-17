#pragma once

#include "iostream" 
#include "cmath"
#include "pch.h"
#include "iomanip"
#include "status.h"

/*using namespace std;*/

	//�������� � ���������
std::vector <double> multAndSum(std::vector<double> xi, std::vector<double> k, double C) {
	std::vector <double> b(k.size());
	for (unsigned int i = 0; i <= k.size() - 1; i++) {
		b[i] = xi[i] + k[i] * C;
	}
	return(b);
}

void runge(status &sv, double h) {
	static std::vector <double> xi;
	xi = sv.getParam();
	//int size = xi.size();
	sv.setParam(xi);
	sv.nonIntegr();

	static std::vector <double> k1;
	k1 = sv.rightPart();
	//static std::vector <double> a(xi.size());
	static std::vector <double> null(xi.size());     //������� ������, �.�. ��� ��������=0
	k1 = k1 * h;
	sv.setParam(k1 * 0.5 + xi);

	sv.nonIntegr();
	static std::vector <double> k2;
	k2 = sv.rightPart();
	k2 = k2 * h;
	sv.setParam(k2 * 0.5 + xi);
	sv.nonIntegr();

	static std::vector <double> k3;
	k3 = sv.rightPart();
	k3 = k3 * h;
	sv.setParam(k3 * 1 + xi);

	sv.nonIntegr();
	static std::vector <double> k4;
	k4 = sv.rightPart();
	k4 = k4 * h;
	static std::vector <double> result(xi.size());

	result = xi + (k1 + (k2 * 2) + (k3 * 2) + k4) * (1.0 / 6.0);
	
	
	sv.setParam(result);
	sv.nonIntegr();
}