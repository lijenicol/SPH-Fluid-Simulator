#pragma once

#include "Core.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

class Geometry
{
private:
	//geometry data
	std::vector<glm::vec3> points;
	std::vector<glm::vec3> normals;
	std::vector<glm::uvec3> triangles;
	std::vector<glm::uvec3> normalIndices;

public:
	Geometry(std::string filename);
	~Geometry();

	void draw(const glm::mat4& modelMtx, const glm::mat4& viewProjMtx, uint shader);
	void init(std::string filename);

	std::vector<unsigned int> indices;
	GLuint vao, vbo, vbo_n, ebo;
};

