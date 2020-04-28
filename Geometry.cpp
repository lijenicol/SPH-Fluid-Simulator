#include "Geometry.h"

const float SCALE_FACTOR = 1.0f;

Geometry::Geometry(std::string filename) {
	init(filename);
}

void Geometry::init(std::string filename) {
	//
	//READ IN GEOMETRY DATA FROM FILE
	//
	std::vector<glm::vec3> inputNormals;
	std::vector<glm::vec3> inputPoints;
	std::ifstream objStream(filename);
	std::string line;
	if (objStream.is_open()) {
		//go through every line of the obj file
		while (getline(objStream, line)) {
			std::stringstream ss;
			ss << line;

			//check for vertex label
			std::string label;
			ss >> label;
			if (label == "v") {
				//construct a point from the line
				glm::vec3 point;
				ss >> point.x;
				ss >> point.y;
				ss >> point.z;

				//push the point to the point list
				inputPoints.push_back(point);
			}
			if (label == "vn") {
				glm::vec3 normal;
				ss >> normal.x;
				ss >> normal.y;
				ss >> normal.z;

				//push the normal to the normals list
				inputNormals.push_back(normal);
			}
			if (label == "f") {
				glm::uvec3 triangle;
				glm::uvec3 normalIndex;
				std::string s;
				int pos;
				std::string v;

				//search face for / delimeter
				ss >> s;
				pos = s.find("/");
				v = s.substr(0, pos);
				triangle.x = std::stoi(v) - 1;
				s.erase(0, pos + 1);
				pos = s.find("/");
				s.erase(0, pos + 1);
				normalIndex.x = std::stoi(s) - 1;

				ss >> s;
				pos = s.find("/");
				v = s.substr(0, pos);
				triangle.y = std::stoi(v) - 1;
				s.erase(0, pos + 1);
				pos = s.find("/");
				s.erase(0, pos + 1);
				normalIndex.y = std::stoi(s) - 1;

				ss >> s;
				pos = s.find("/");
				v = s.substr(0, pos);
				triangle.z = std::stoi(v) - 1;
				s.erase(0, pos + 1);
				pos = s.find("/");
				s.erase(0, pos + 1);
				normalIndex.z = std::stoi(s) - 1;

				triangles.push_back(triangle);
				normalIndices.push_back(normalIndex);
			}
		}
		objStream.close();
	}
	else {
		std::cout << "File cannot be read or does not exist: " << filename;
	}

	//reorder and duplicate indices so that indices line up (thanks normals)
	for (unsigned int i = 0; i < triangles.size(); i++) {
		points.push_back(inputPoints[triangles[i].x]);
		points.push_back(inputPoints[triangles[i].y]);
		points.push_back(inputPoints[triangles[i].z]);
		normals.push_back(inputNormals[normalIndices[i].x]);
		normals.push_back(inputNormals[normalIndices[i].y]);
		normals.push_back(inputNormals[normalIndices[i].z]);
		indices.push_back(i * 3);
		indices.push_back(i * 3 + 1);
		indices.push_back(i * 3 + 2);
	}

	//
	//SETUP OPENGL BUFFERS HERE
	//

	// Generate a vertex array (VAO) and a vertex buffer objects (VBO) + EBO.
	glGenVertexArrays(1, &vao);
	glGenBuffers(1, &vbo);
	glGenBuffers(1, &vbo_n);
	glGenBuffers(1, &ebo);

	// Bind to the VAO.
	glBindVertexArray(vao);

	// Bind to the first VBO. We will use it to store the points.
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	// Pass in the data.
	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * points.size(),
		points.data(), GL_STATIC_DRAW);
	// Enable vertex attribute 0. 
	// We will be able to access points through it.
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

	//Send normal data to gpu
	glBindBuffer(GL_ARRAY_BUFFER, vbo_n);
	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * normals.size(), normals.data(), GL_STATIC_DRAW);

	//Bind to the EBO 
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	//Send the information to the gpu
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * indices.size(), indices.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

	// Unbind from the VBO.
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	// Unbind from the VAO.
	glBindVertexArray(0);
}

Geometry::~Geometry() {
	// Delete the VBO and the VAO.
	glDeleteBuffers(1, &vbo);
	glDeleteBuffers(1, &vbo_n);
	glDeleteBuffers(1, &ebo);
	glDeleteVertexArrays(1, &vao);
}

void Geometry::draw(const glm::mat4& modelMtx, const glm::mat4& viewProjMtx, uint shader)
{
	glUseProgram(shader);

	// Set model view matrix
	glUniformMatrix4fv(glGetUniformLocation(shader, "ModelMtx"), 1, false, (float*)&modelMtx);
	glm::mat4 mvpMtx = viewProjMtx * modelMtx;
	glUniformMatrix4fv(glGetUniformLocation(shader, "ModelViewProjMtx"), 1, false, (float*)&mvpMtx);

	// Bind to the VAO.
	glBindVertexArray(vao);
	// Draw the model
	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
	// Unbind from the VAO.
	glBindVertexArray(0);

	glUseProgram(0);
}