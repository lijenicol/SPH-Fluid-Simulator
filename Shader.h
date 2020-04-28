////////////////////////////////////////
// Shader.h
////////////////////////////////////////

#pragma once

#include "Core.h"

////////////////////////////////////////////////////////////////////////////////

class Shader {
public:
	enum ShaderType {eGeometry,eVertex,eFragment,eCompute};

	Shader(const char *filename,ShaderType type);
	~Shader();

	// Access functions
	uint GetShaderID()				{return ShaderID;}

private:
	uint ShaderID;
};

////////////////////////////////////////////////////////////////////////////////

class ShaderProgram {
public:
	enum ProgramType {eGeometry,eRender,eCompute};

	ShaderProgram(const char *filename,ProgramType type);
	~ShaderProgram();

	// Access functions
	uint GetProgramID() const				{return ProgramID;}

private:
	enum {eMaxShaders=3};
	ProgramType Type;
	Shader *Shaders[eMaxShaders];
	uint ProgramID;
};

////////////////////////////////////////////////////////////////////////////////
