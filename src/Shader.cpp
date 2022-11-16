////////////////////////////////////////
// Shader.cpp
////////////////////////////////////////

#include "Shader.h"

#include <string>
#include <vector>
#include <fstream>

////////////////////////////////////////////////////////////////////////////////
// Shader
////////////////////////////////////////////////////////////////////////////////

Shader::Shader(const char *filename,ShaderType type) {
	// Open shader file
	std::string fileData = "";
	std::fstream stream(filename, std::ios::in);
	if(!stream.is_open()) {
		printf("ERROR: Can't open shader file '%s'\n",filename);
		ShaderID=0;
		return;
	}

	// Read shader file
	std::string line = "";
	while (getline(stream, line)) {
		fileData += "\n" + line;
	}

	// Close file
	stream.close();

	// Append header & create shader
	std::string shaderSource;
	switch(type) {
		case eGeometry:
			shaderSource="#define GEOMETRY_SHADER\n"+fileData;
			ShaderID=glCreateShader(GL_GEOMETRY_SHADER);
			break;
		case eVertex:
			shaderSource="#define VERTEX_SHADER\n"+fileData;
			ShaderID=glCreateShader(GL_VERTEX_SHADER);
			break;
		case eFragment:
			shaderSource="#define FRAGMENT_SHADER\n"+fileData;
			ShaderID=glCreateShader(GL_FRAGMENT_SHADER);
			break;
		case eCompute:
			shaderSource="#define COMPUTE_SHADER\n"+fileData;
			ShaderID=glCreateShader(GL_COMPUTE_SHADER);
			break;
	}
	shaderSource="#version 430\n"+shaderSource;

	// Compile shader
	const char* rawShaderSource = shaderSource.c_str();
	glShaderSource(ShaderID, 1, &rawShaderSource, NULL);
	glCompileShader(ShaderID);

	// Verify compile worked
	GLint isCompiled = 0;
	glGetShaderiv(ShaderID, GL_COMPILE_STATUS, &isCompiled);
	if (isCompiled == GL_FALSE) {
		printf("ERROR: Can't compile shader '%s'\n",filename);

		// Print error message
		GLint maxLength = 0;
		glGetShaderiv(ShaderID, GL_INFO_LOG_LENGTH, &maxLength);
		std::vector<GLchar> errorLog(maxLength);
		glGetShaderInfoLog(ShaderID, maxLength, &maxLength, &errorLog[0]);
		printf("GL ERRORS:\n%s\n",&errorLog[0]);

		// Delete shader
		glDeleteShader(ShaderID);
		ShaderID=0;
	}
}

////////////////////////////////////////////////////////////////////////////////

Shader::~Shader() {
	glDeleteShader(ShaderID);
}

////////////////////////////////////////////////////////////////////////////////
// ShaderProgram
////////////////////////////////////////////////////////////////////////////////

ShaderProgram::ShaderProgram(const char *filename,ProgramType type) {
	Type=type;
	for(int i=0;i<eMaxShaders;i++)
		Shaders[i]=0;

	// Load & compile shaders
	if(Type==eGeometry) {
		Shaders[0]=new Shader(filename,Shader::eVertex);
		Shaders[1]=new Shader(filename,Shader::eGeometry);
		Shaders[2]=new Shader(filename,Shader::eFragment);
	}
	else if(Type==eRender) {
		Shaders[0]=new Shader(filename,Shader::eVertex);
		Shaders[1]=new Shader(filename,Shader::eFragment);
	}
	else if(Type==eCompute) {
		Shaders[0]=new Shader(filename,Shader::eCompute);
	}

	// Link shader program
	ProgramID=glCreateProgram();
	for(int i=0;i<eMaxShaders;i++) {
		if(Shaders[i]==0)
			break;
		glAttachShader(ProgramID,Shaders[i]->GetShaderID());
	}
	glLinkProgram(ProgramID);
}

////////////////////////////////////////////////////////////////////////////////

ShaderProgram::~ShaderProgram() {
	for(int i=0;i<eMaxShaders;i++)
		delete Shaders[i];
	glDeleteProgram(ProgramID);
}

////////////////////////////////////////////////////////////////////////////////
