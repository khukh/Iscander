#pragma once
class targetStatus
{
	double X, Y, Z;
	double Vx, Vy, Vz;
public:
	targetStatus(std::vector<double> a);
	~targetStatus();
};

