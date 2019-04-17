#include "pch.h"
#include "rotation.h"





rotation::rotation(double pitch, double yaw, double roll):A(pitch,yaw,roll),RG(pitch, yaw, roll) {
	Angles.resize(3);
	Angles[0] = pitch;
	Angles[0] = yaw;
	Angles[0] = roll;
}

rotation::~rotation() {}
