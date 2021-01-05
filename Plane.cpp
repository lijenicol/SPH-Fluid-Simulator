#include "Plane.h"

Plane::Plane() {
	//vertices of plane
	float vertices[] = {
	 0.5f,  0.f, 0.5f,  // top right
	 0.5f, -0.f, -0.5f,  // bottom right
	-0.5f, -0.f, -0.5f,  // bottom left
	-0.5f,  0.f, 0.5f   // top left 
	};

	//extend vertices of plane by some multiple
	for (int i = 0; i < 12; i++) {
		vertices[i] *= 1000.0f;
	}

	//normals of plane
	float normals[] = {
	 0.f,  1.f, 0.f,  // top right
	 0.f,  1.f, -0.f,  // bottom right
	-0.f,  1.f, -0.f,  // bottom left
	-0.f,  1.f, 0.f   // top left 
	};

	//indices of plane
	unsigned int indices[] = { 
		0, 1, 3,   // first triangle
		1, 2, 3    // second triangle
	};

	//generate buffers
	glGenVertexArrays(1, &vao);
	glGenBuffers(1, &vbo);
	glGenBuffers(1, &vbon);
	glGenBuffers(1, &ebo);

	// Bind to the VAO.
	glBindVertexArray(vao);

	// Send vertex data
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glEnableVertexAttribArray(0);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices),
		vertices, GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

	//Send normal data to gpu
	glBindBuffer(GL_ARRAY_BUFFER, vbon);
	glEnableVertexAttribArray(1);
	glBufferData(GL_ARRAY_BUFFER, sizeof(normals), normals, GL_STATIC_DRAW);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

	//Send ebo data 
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	// Unbind from the VBO.
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	// Unbind from the VAO.
	glBindVertexArray(0);
}

void Plane::draw(glm::mat4 viewProjMtx, GLuint shader) {
	glm::mat4 modelMtx = glm::mat4(1);

	// Set up shader
	glUseProgram(shader);
	glUniformMatrix4fv(glGetUniformLocation(shader, "ModelMtx"), 1, false, (float*)&modelMtx);

	glm::mat4 mvpMtx = viewProjMtx;
	glUniformMatrix4fv(glGetUniformLocation(shader, "viewProjMtx"), 1, false, (float*)&mvpMtx);

	// Bind to the VAO.
	glBindVertexArray(vao);
	// Draw the model
	glDisable(GL_CULL_FACE);	//disable face culling so that both sides of plane are rendered
	glDrawElements(GL_TRIANGLE_STRIP, 6, GL_UNSIGNED_INT, 0);
	glEnable(GL_CULL_FACE);
	// Unbind from the VAO.
	glBindVertexArray(0);

	glUseProgram(0);
}