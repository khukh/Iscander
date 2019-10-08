#pragma once
const double PI = 3.14159265358979323846264338328;
const double toDeg = 180 / PI;
const double toRad = PI / 180;

const int N = 6;

const double X0 = 0;
const double Y0 = 10000 - 150 * N;
const double Z0 = 0;

const double ALPHA0 = 0;

const double PITCH0 = (-60 + 0.5 * N) * toRad; //������ � ��������
const double ROLL0 = (0)*toRad; //���� ����� � ��������
const double YAW0 = (0)*toRad; //���� �������� � ��������

const double W_X0 = 0;
const double W_Y0 = 0;
const double W_Z0 = 0;

const double I_X0 = 170;
const double I_Y0 = 640;
const double I_Z0 = 640;

const double D_M = 0.95;
const double L = 7;

const double V0 = 1200; //������ ��������
const double M0 = 1500 - 5 * N;

const double H = 0.01;
const double T_FIN = 10.0;
const double EPS1 = 1E-6;


const double K_T1 = 1.0;
const double K_N1 = 1.0;
const double K_E1 = 1.0;
const double K_E2 = 1.0;


