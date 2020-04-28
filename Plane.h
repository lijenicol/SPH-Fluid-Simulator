#pragma once
#include "GL/glew.h"
#include <glm/glm.hpp>
#include <vector>

class Plane
{
private:
	//geometry properties
	GLuint vao, vbo, vbon, ebo;
public:
	Plane();
	~Plane();

	void draw(glm::mat4 viewProjMtx, GLuint shader);
};

