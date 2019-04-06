#pragma once

#include "iostream" 
#include "cmath"
#include "vector"
#include "iomanip"
#include "status.h"

/*using namespace std;*/

	//действия с векторами
std::vector <double> multAndSum(std::vector<double> xi, std::vector<double> k, double C) {
	std::vector <double> b(k.size());
	for (unsigned int i = 0; i <= k.size() - 1; i++) {
		b[i] = xi[i] + k[i] * C;
	}
	return(b);
}

void runge(status &sv, double h, int n) {
	static std::vector <double> xi;
	xi = sv.getParam();
	//int size = xi.size();
	sv.setParam(xi);
	sv.nonIntegr( n);

	static std::vector <double> k1;
	k1 = sv.rightPart(n);
	//static std::vector <double> a(xi.size());
	static std::vector <double> null(xi.size());     //нулевой вектор, т.е. все элементы=0
	k1 = multAndSum(null, k1, h);
	sv.setParam(multAndSum(xi, k1, 0.5));

	sv.nonIntegr( n);
	static std::vector <double> k2;
	k2 = sv.rightPart(n);
	k2 = multAndSum(null, k2, h);
	sv.setParam(multAndSum(xi, k2, 0.5));
	sv.nonIntegr(n);

	static std::vector <double> k3;
	k3 = sv.rightPart(n);
	k3 = multAndSum(null, k3, h);
	sv.setParam(multAndSum(xi, k3, 1));

	sv.nonIntegr(n);
	static std::vector <double> k4;
	k4 = sv.rightPart(n);
	k4 = multAndSum(null, k4, h);
	static std::vector <double> result(xi.size());
	
	for (unsigned int i = 0; i <= xi.size() - 1; i++) {
		result[i] = xi[i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
	}
	sv.setParam(result);
	sv.nonIntegr(n);
}