////////////////////////////////////////
// Tester.cpp
////////////////////////////////////////

#include "Tester.h"
#include <iostream>

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glut.h"
#include "imgui/imgui_impl_opengl2.h"

////////////////////////////////////////////////////////////////////////////////

static Tester *TESTER=0;

int main(int argc, char **argv) {
	glutInit(&argc, argv);

	TESTER=new Tester("SPH Fluid Simulator",argc,argv);
	glutMainLoop();
	delete TESTER;

	return 0;
}

////////////////////////////////////////////////////////////////////////////////

// These are really HACKS to make glut call member functions instead of static functions
static void display()									{TESTER->Draw();}
static void idle()										{TESTER->Update();}
static void resize(int x,int y)							{TESTER->Resize(x,y);}
static void keyboard(unsigned char key,int x,int y)		{TESTER->Keyboard(key,x,y);}
static void specialKeys(int key, int x, int y) { TESTER->SpecialKeys(key, x, y); }
static void mousebutton(int btn,int state,int x,int y)	{TESTER->MouseButton(btn,state,x,y);}
static void mousemotion(int x, int y)					{TESTER->MouseMotion(x,y);}

////////////////////////////////////////////////////////////////////////////////

Tester::Tester(const char *windowTitle,int argc,char **argv) {
	WinX=800;
	WinY=600;
	LeftDown=MiddleDown=RightDown=false;
	MouseX=MouseY=0;
	prevTime = 0;
	currentTime = 0;
	deltaTime = 0;

	// Create the window
	glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );
	glutInitWindowSize( WinX, WinY );
	glutInitWindowPosition( 100, 100 );
	WindowHandle = glutCreateWindow( windowTitle );
	glutSetWindowTitle( windowTitle );
	glutSetWindow( WindowHandle );

	// Background color
	glClearColor( 0., 0., 0., 1. );

	// Callbacks
	glutDisplayFunc( display );
	glutIdleFunc( idle );
	glutKeyboardFunc( keyboard );
	glutMouseFunc( mousebutton );
	glutMotionFunc( mousemotion );
	glutPassiveMotionFunc( mousemotion );
	glutReshapeFunc( resize );

	// Initialize GLEW
	glewInit();
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	// Setup Dear ImGui context
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;

	// Setup Dear ImGui style
	ImGui::StyleColorsDark();

	// Setup Platform/Renderer bindings
	ImGui_ImplGLUT_Init();
	ImGui_ImplGLUT_InstallFuncs();
	ImGui_ImplOpenGL2_Init();

	// Initialize components
	Program=new ShaderProgram("Model.glsl",ShaderProgram::eRender);
	Cam=new Camera;
	Cam->SetAspect(float(WinX)/float(WinY));

	//initialize plane (only for rendering)
	plane = new Plane();

	//init SPH system
	sphSystem = new SPHSystem(5, 0.02f, 1000, 0.04f, 0.04f, -9.8f, 0.2f);
}

////////////////////////////////////////////////////////////////////////////////

Tester::~Tester() {
	delete Program;
	delete Cam;

	// Cleanup imgui
	ImGui_ImplOpenGL2_Shutdown();
	ImGui_ImplGLUT_Shutdown();
	ImGui::DestroyContext();

	glFinish();
	glutDestroyWindow(WindowHandle);
}

////////////////////////////////////////////////////////////////////////////////

void Tester::Update() {
	//calculate delta time
	currentTime = glutGet(GLUT_ELAPSED_TIME);
	deltaTime = (currentTime - prevTime)/1000.f;
	prevTime = currentTime;

	// Update the components in the world
	Cam->Update();
	
	//update sph system
	sphSystem->update(deltaTime);

	// Tell glut to re-display the scene
	glutSetWindow(WindowHandle);
	glutPostRedisplay();
}

////////////////////////////////////////////////////////////////////////////////

void Tester::Reset() {
	Cam->Reset();
	Cam->SetAspect(float(WinX)/float(WinY));
}

////////////////////////////////////////////////////////////////////////////////

void Tester::Draw() {
	// Begin drawing scene
	glViewport(0, 0, WinX, WinY);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Draw plane
	plane->draw(Cam->GetViewProjectMtx(), Program->GetProgramID());

	// Draw sph system
	sphSystem->draw(Cam->GetViewProjectMtx(), Program->GetProgramID());

	// Render GUI
	ImGui_ImplOpenGL2_NewFrame();
	ImGui_ImplGLUT_NewFrame();
	{
		static int numParticles = 5;
		static float nMass = 0.065f;
		static float nh = 0.065f;
		static float nRest = 1000.f;
		static float nVisco = 3.5f;
		static int counter = 0;

		ImGui::Begin("SPH debug");                          // Create GUI window

		ImGui::Text("Change values for the simulation. Press RESET to commit changes"); 

		ImGui::SliderInt("Number of Particles", &numParticles, 10, 600);            // Edit number of particles
		ImGui::SliderFloat("Mass of Particles", &nMass, 0.001f, 1.f);            // Edit mass
		ImGui::SliderFloat("Support Radius", &nh, 0.001f, 1.f);            // Edit support radius
		ImGui::SliderFloat("Rest Density", &nRest, 0.001f, 2000.f);            // Edit rest density
		ImGui::SliderFloat("Viscosity Constant", &nVisco, 0.001f, 5.f);            // Edit viscosity

		if (ImGui::Button("RESET")) {
			delete sphSystem;
			sphSystem = new SPHSystem(numParticles, nMass, nRest, nVisco, nh, -9.8, 1.f);
		}

		ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
		ImGui::End();
	}
	glUseProgram(0);
	ImGui::Render();
	ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());

	// Finish drawing scene
	glFinish();
	glutSwapBuffers();
}

////////////////////////////////////////////////////////////////////////////////

void Tester::Quit() {
	glFinish();
	glutDestroyWindow(WindowHandle);
	exit(0);
}

////////////////////////////////////////////////////////////////////////////////

void Tester::Resize(int x,int y) {
	WinX = x;
	WinY = y;
	Cam->SetAspect(float(WinX)/float(WinY));

	ImGui_ImplGLUT_ReshapeFunc(x, y);
}

////////////////////////////////////////////////////////////////////////////////

void Tester::Keyboard(int key,int x,int y) {
	ImGui_ImplGLUT_KeyboardFunc(key, x, y);

	switch(key) {
		case 0x1b:		// Escape
			Quit();
			break;
		case 'r':
			Reset();
			sphSystem->reset();
			break;
		case 'p':
			sphSystem->print();
			break;
		case 'c':
			sphSystem->startSimulation();
			break;
	}
}

////////////////////////////////////////////////////////////////////////////////

void Tester::SpecialKeys(int key, int x, int y) {
	switch (key) {	
		
	}
}

////////////////////////////////////////////////////////////////////////////////

void Tester::MouseButton(int btn,int state,int x,int y) {
	// Send mouse inputs to GUI
	ImGui_ImplGLUT_MouseFunc(btn, state, x, y);

	// Don't read input if mouse inside GUI
	if (ImGui::GetIO().WantCaptureMouse)
		return;

	if(btn==GLUT_LEFT_BUTTON) {
		LeftDown = (state==GLUT_DOWN);
	}
	else if(btn==GLUT_MIDDLE_BUTTON) {
		MiddleDown = (state==GLUT_DOWN);
	}
	else if(btn==GLUT_RIGHT_BUTTON) {
		RightDown = (state==GLUT_DOWN);
	}
}

////////////////////////////////////////////////////////////////////////////////

void Tester::MouseMotion(int nx,int ny) {
	// Send mouse inputs to GUI
	ImGui_ImplGLUT_MotionFunc(nx, ny);

	// Don't read input if mouse inside GUI
	if (ImGui::GetIO().WantCaptureMouse)
		return;

	int maxDelta=100;
	int dx = glm::clamp(nx - MouseX,-maxDelta,maxDelta);
	int dy = glm::clamp(-(ny - MouseY),-maxDelta,maxDelta);

	MouseX = nx;
	MouseY = ny;

	// Move camera
	// NOTE: this should really be part of Camera::Update()
	if(LeftDown) {
		const float rate=1.0f;
		Cam->SetAzimuth(Cam->GetAzimuth()+dx*rate);
		Cam->SetIncline(glm::clamp(Cam->GetIncline()-dy*rate,-90.0f,90.0f));
	}
	if(RightDown) {
		const float rate=0.005f;
		float dist=glm::clamp(Cam->GetDistance()*(1.0f-dx*rate),0.01f,1000.0f);
		Cam->SetDistance(dist);
	}
}

////////////////////////////////////////////////////////////////////////////////
