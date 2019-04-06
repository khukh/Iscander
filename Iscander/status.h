#pragma once

#include "pch.h"

#include "iostream" 
#include <vector>
#include "Constants.h"
#include <fstream>

#include <RGParam.h>



class status {
public:
	status();
	~status();

	RGParam Rg;



	void nonIntegr(int n);
	void printParam( std::ofstream &fout);	//����� ����������
	void setParam(std::vector <double> b);
	void setPSI(std::vector <std::vector <double>> psiX, std::vector <std::vector <double>> psiV);
	void setUcontrol(int n);
	std::vector <double> getParam();
	std::vector <double> rightPart(int n);	//�������� �����������	
	double r;// ������-������
	status& operator=(const status& right);
protected:
	std::vector<double> parametr;
	/*	parametr [0] = Vx
		parametr [1] = Vy
		parametr [2] = Vz
		parametr [3] = x
		parametr [4] = y
		parametr [5] = z
		parametr [6] = t
		parametr [7] = �������� (Ux)^2
		parametr [8] = �������� (Uy)^2
		parametr [9] = �������� (Uz)^2
		*/
	std::vector <std::vector <double>> psiX;
	std::vector <std::vector <double>> psiV;
	std::vector<double> u;
	
	double mu;


};

