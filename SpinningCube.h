////////////////////////////////////////
// SpinningCube.h
////////////////////////////////////////

#pragma once

#include "Model.h"

////////////////////////////////////////////////////////////////////////////////

// SpinningCube is an example of a basic animating object. It can be used as a
// pattern for creating more complex objects.

class SpinningCube {
public:
	SpinningCube();

	void Update();
	void Reset();
	void Draw(const glm::mat4 &viewProjMtx,uint shader);

private:
	// Constants
	Model CubeModel;
	glm::vec3 Position;
	glm::vec3 Axis;
	float SpinDelta;

	// Variables
	float Angle;
	glm::mat4x4 WorldMtx;
};

////////////////////////////////////////////////////////////////////////////////
