#include "pch.h"
#include "KoefCalc.h"


/////////////////////
double CyDeltaT(double mach, double alpha) {
	double alDeg = alpha * toDeg;
	double p1 = 1 / (243.84E-3 / exp(alDeg) + 74.309E-3);
	double p2 = log10(1.9773*alDeg*alDeg - 25.587*alDeg + 83.354);
	double p3 = 18.985*alDeg*alDeg - 375.76*alDeg + 1471;
	double p4 = -51.164E-3*alDeg*alDeg + 805.52E-3*alDeg + 1.8929;
	double k = 2 * (-p1 * 1E-6*mach*mach + p2 * 1E-12*exp(mach) - p3 * 1E-6*mach - p4 * 1E-3);
	return k;
}

double CyAl(double mach, double alpha) {
	double ds = 1.86 * (11.554 / exp(mach) - (2.5191E-3)*mach*mach - 5.024 / mach + (52.836E-3)*mach + 4.112);
	if (ds >= 0) {
		return sqrt(ds);
	}
	else {
		return 1.86*1.039;
	}
}

double CzDeltaN(double mach, double betta) {
	return -CyDeltaT(mach, betta);
}

double CzBetta(double mach, double betta) {
	return -CyAl(mach, betta);
}

///////////////////////

double Cx(double mach, double alphaSpace) {
	double cx = 1 / (73.211 / exp(mach) - 47.483 / mach + 16.878);
	return cx;
}

double Cy(double mach, double alpha, double deltaT) {
	double cy = CyAl(mach, alpha) * alpha + CyDeltaT(mach, alpha) * deltaT;
	return cy;
}
double Cz(double mach, double betta, double deltaN) {
	double cz = CzBetta(mach, betta) * betta + CzDeltaN(mach, betta) * deltaN;
	return cz;
}

////////////////////////
double MxOmegaX(double mach, double alpha) {
	return -0.005*0.6786;

}

double MxDelta(double mach, double alpha, double q) {
	double p = 1000;
	return -p / q;
}
/////////
double MzOmegaZ(double mach, double alpha) {
	double m = 1.89*(146.79E-6*mach*mach - 158.98E-3 / mach - 7.6396E-3*mach - 68.195E-3);
	return m;

}

double MzAlpha(double mach, double alpha) {
	double m = (-766.79E-3/exp(mach) + 438.74E-3 / mach + 5.8822E-3*mach - 158.34E-3);
	return m;

}

double MzDelta(double mach, double alpha, double q) {
	double k1 = 0;
	double k2 = 0;
	double a = abs(alpha*toDeg);
	k1 = exp(-19.488E-3*a*a - 378.62E-3*a + 6.7518);
	k2 = exp(-21.234E-3*a*a - 635.84E-6*exp(a) - 98.296E-3*a + 2.5938);//?????????????????????
	double res = 1.89*sqrt(k1*1E-9*mach*mach+k2*1E-6);
	return res;
}
/////////
double MyOmegaY(double mach, double betta) {
	return MzOmegaZ(mach, betta);

}

double MyBetta(double mach, double betta) {
	return MzAlpha(mach, betta);
}

double MyDelta(double mach, double betta, double q) {
	return MzDelta(mach, betta, q);
}