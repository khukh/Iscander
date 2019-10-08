#pragma once
double Cx(double mach, double alphaSpace);
double Cy(double mach, double alpha, double deltaT);
double Cz(double mach, double betta, double deltaN);

double CyDeltaT(double mach, double alpha);
double CyAl(double mach, double alpha);

double CzDeltaN(double mach, double betta);
double CzBetta(double mach, double betta);

////////////////////////
double MxOmegaX(double mach, double alpha);
double MxDelta(double mach, double alpha, double q);

double MzOmegaZ(double mach, double alpha);
double MzAlpha(double mach, double alpha);
double MzDelta(double mach, double alpha, double q);

double MyOmegaY(double mach, double betta);
double MyBetta(double mach, double betta);
double MyDelta(double mach, double betta, double q);