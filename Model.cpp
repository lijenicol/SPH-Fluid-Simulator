////////////////////////////////////////
// Model.cpp
////////////////////////////////////////

#include "Model.h"

////////////////////////////////////////////////////////////////////////////////

Model::Model() {
	glGenBuffers(1,&VertexBuffer);
	glGenBuffers(1,&IndexBuffer);

	Count=0;
}

////////////////////////////////////////////////////////////////////////////////

Model::~Model() {
	glDeleteBuffers(1,&IndexBuffer);
	glDeleteBuffers(1,&VertexBuffer);
}

////////////////////////////////////////////////////////////////////////////////

void Model::Draw(const glm::mat4 &modelMtx,const glm::mat4 &viewProjMtx,uint shader) {
	// Set up shader
	glUseProgram(shader);
	glUniformMatrix4fv(glGetUniformLocation(shader,"ModelMtx"),1,false,(float*)&modelMtx);

	glm::mat4 mvpMtx=viewProjMtx * modelMtx;
	glUniformMatrix4fv(glGetUniformLocation(shader,"ModelViewProjMtx"),1,false,(float*)&mvpMtx);

	// Set up state
	glBindBuffer(GL_ARRAY_BUFFER,VertexBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,IndexBuffer);

	uint posLoc=0;
	glEnableVertexAttribArray(posLoc);
	glVertexAttribPointer(posLoc,3,GL_FLOAT,GL_FALSE,sizeof(ModelVertex),0);

	uint normLoc=1;
	glEnableVertexAttribArray(normLoc);
	glVertexAttribPointer(normLoc,3,GL_FLOAT,GL_FALSE,sizeof(ModelVertex),(void*)12);

	// Draw geometry
	glDrawElements(GL_TRIANGLES,Count,GL_UNSIGNED_INT,0);

	// Clean up state
	glDisableVertexAttribArray(normLoc);
	glDisableVertexAttribArray(posLoc);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
	glBindBuffer(GL_ARRAY_BUFFER,0);

	glUseProgram(0);
}

////////////////////////////////////////////////////////////////////////////////

void Model::MakeBox(const glm::vec3 &boxMin,const glm::vec3 &boxMax) {
	// Specify vertices
	std::vector<ModelVertex> vtx={
		// Front
		{{boxMin.x,boxMin.y,boxMax.z},	{0,0,1}},
		{{boxMax.x,boxMin.y,boxMax.z},	{0,0,1}},
		{{boxMax.x,boxMax.y,boxMax.z},	{0,0,1}},
		{{boxMin.x,boxMax.y,boxMax.z},	{0,0,1}},

		// Back
		{{boxMax.x,boxMin.y,boxMin.z},	{0,0,-1}},
		{{boxMin.x,boxMin.y,boxMin.z},	{0,0,-1}},
		{{boxMin.x,boxMax.y,boxMin.z},	{0,0,-1}},
		{{boxMax.x,boxMax.y,boxMin.z},	{0,0,-1}},

		// Top
		{{boxMin.x,boxMax.y,boxMax.z},	{0,1,0}},
		{{boxMax.x,boxMax.y,boxMax.z},	{0,1,0}},
		{{boxMax.x,boxMax.y,boxMin.z},	{0,1,0}},
		{{boxMin.x,boxMax.y,boxMin.z},	{0,1,0}},

		// Bottom
		{{boxMin.x,boxMin.y,boxMin.z},	{0,-1,0}},
		{{boxMax.x,boxMin.y,boxMin.z},	{0,-1,0}},
		{{boxMax.x,boxMin.y,boxMax.z},	{0,-1,0}},
		{{boxMin.x,boxMin.y,boxMax.z},	{0,-1,0}},

		// Left
		{{boxMin.x,boxMin.y,boxMin.z},	{-1,0,0}},
		{{boxMin.x,boxMin.y,boxMax.z},	{-1,0,0}},
		{{boxMin.x,boxMax.y,boxMax.z},	{-1,0,0}},
		{{boxMin.x,boxMax.y,boxMin.z},	{-1,0,0}},

		// Right
		{{boxMax.x,boxMin.y,boxMax.z},	{1,0,0}},
		{{boxMax.x,boxMin.y,boxMin.z},	{1,0,0}},
		{{boxMax.x,boxMax.y,boxMin.z},	{1,0,0}},
		{{boxMax.x,boxMax.y,boxMax.z},	{1,0,0}},
	};

	// Specify indices
	std::vector<uint> idx={
		0,1,2,		0,2,3,			// Front
		4,5,6,		4,6,7,			// Back
		8,9,10,		8,10,11,		// Top
		12,13,14,	12,14,15,		// Bottom
		16,17,18,	16,18,19,		// Left
		20,21,22,	20,22,23,		// Right
	};

	// Create vertex & index buffers
	SetBuffers(vtx,idx);
}

////////////////////////////////////////////////////////////////////////////////

void Model::SetBuffers(const std::vector<ModelVertex> &vtx,const std::vector<uint> &idx) {
	Count=(int)idx.size();

	// Store vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER,VertexBuffer);
	glBufferData(GL_ARRAY_BUFFER,vtx.size()*sizeof(ModelVertex),&vtx[0],GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER,0);

	// Store index buffer
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,IndexBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,idx.size()*sizeof(uint),&idx[0],GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
}

////////////////////////////////////////////////////////////////////////////////
