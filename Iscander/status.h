#pragma once

#include "pch.h"

#include "iostream" 
//#include <vector>
//#include "Constants.h"
#include <fstream>
#include "matrix.h"
#include "RGParam.h"



class status {
public:
	status();
	~status();

	

	std::vector <double> ForcePr;
	std::vector <double> ForcePrG;  //�������� 
	std::vector <double> ForceG;  //������ ���� ����������
	std::vector <double> Torque;  //������ ����

	matrix A;
	RGParam Rg;

	void nonIntegr();  //�������� ��������������� ����������
	std::vector <double> rightPart();	//�������� �����������	

	void printParam( std::ofstream &fout);	//����� ����������
	void setParam(std::vector <double> b);
	std::vector <double> getParam();


	status& operator=(const status& right);
protected:
	std::vector<double> parametr;
	/*	parametr [0] = Vx
		parametr [1] = Vy
		parametr [2] = Vz
		parametr [3] = x
		parametr [4] = y
		parametr [5] = z
		parametr [6] = Wx
		parametr [7] = Wy
		parametr [8] = Wz
		parametr [9] = ro
		parametr [9] = ro
		parametr [9] = ro
		parametr [9] = ro

		*/
	//std::vector <std::vector <double>> psiX;
	//std::vector <std::vector <double>> psiV;
	//std::vector<double> u;
	double alpha;
	double betta;
	const double S_M = 0.1;
	const double M = M0;
	const double I_X = I_X0;
	const double I_Y = I_Y0;
	const double I_Z = I_Z0;
	/////TODO: �������� ������ ��� ����������� ����������
	double g = 9.80665;
	double Cx = 0.5;
	double Cy = 0.5;
	double Cz = 0.5;
	double density = 0.8;
	
	//double mu;


};

