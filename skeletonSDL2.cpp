//DH2323 skeleton code, Lab3 (SDL2 version)
#include <iostream>
#include <glm/glm.hpp>
#include "SDL2auxiliary.h"
#include "TestModel.h"
#include <algorithm> //for max()

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::ivec2;

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES
const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL2Aux *sdlAux;
int t;
vector<Triangle> triangles;
vec3 cameraPos(0, 0, -3.001); // Fixed camera position
float focalLength = SCREEN_HEIGHT;
float yaw = 0.0f;   // Horizontal angle (around Y-axis)
float pitch = 0.0f; // Vertical angle (around X-axis)
float roll = 0.0f;
float rotationSpeed = 0.002f; // Radians per ms
float moveSpeed = 0.005f;
// ----------------------------------------------------------------------------
// FUNCTIONS
void Update(void);
void Draw(void);
void VertexShader(const vec3& v, ivec2& p);

int SDL_main(int argc, char* argv[])
{
	LoadTestModel(triangles);  // Load model
	sdlAux = new SDL2Aux(SCREEN_WIDTH, SCREEN_HEIGHT);
	t = SDL_GetTicks();	// Set start value for timer.

	while (!sdlAux->quitEvent())
	{
		Update();
		Draw();
	}
	sdlAux->saveBMP("screenshot.bmp");
	return 0;
}

void Update(void)
{
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

	const Uint8* keystate = SDL_GetKeyboardState(NULL);
	float move = moveSpeed * dt;
	float rot = rotationSpeed * dt;

	// ─── Movement ───
	// Get forward, right, up from yaw and pitch (no roll for movement)
	vec3 forward(sin(yaw), 0, cos(yaw));
	vec3 right(cos(yaw), 0, -sin(yaw));
	vec3 up(0, 1, 0);

	if (keystate[SDL_SCANCODE_UP])    cameraPos += move * forward;
	if (keystate[SDL_SCANCODE_DOWN])  cameraPos -= move * forward;
	if (keystate[SDL_SCANCODE_LEFT])  cameraPos -= move * right;
	if (keystate[SDL_SCANCODE_RIGHT]) cameraPos += move * right;

	// ─── Rotation ───
	if (keystate[SDL_SCANCODE_W]) pitch -= rot;
	if (keystate[SDL_SCANCODE_S]) pitch += rot;
	if (keystate[SDL_SCANCODE_A]) roll -= rot;
	if (keystate[SDL_SCANCODE_D]) roll += rot;
	if (keystate[SDL_SCANCODE_Q]) yaw -= rot;
	if (keystate[SDL_SCANCODE_E]) yaw += rot;

	// Clamp pitch to avoid flipping
	float limit = glm::radians(89.0f);
	pitch = std::max(-limit, std::min(limit, pitch));
}

void Draw()
{
	sdlAux->clearPixels();

	for (int i = 0; i<triangles.size(); ++i)
	{
		vector<vec3> vertices(3);

		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;

		for (int v = 0; v < 3; ++v)
		{
			ivec2 projPos;
			VertexShader(vertices[v], projPos);
			vec3 color(1, 1, 1);
			if (projPos.x >= 0 && projPos.x < SCREEN_WIDTH && projPos.y >= 0 && projPos.y < SCREEN_HEIGHT)
			{
				sdlAux->putPixel(projPos.x, projPos.y, color);
			}
		}
	}

	sdlAux->render();
}

void VertexShader(const vec3& v, ivec2& p)
{
	vec3 relative = v - cameraPos;

	// Build full rotation matrix: R = Rx * Ry * Rz
	mat3 Rx = mat3(1, 0, 0,
				   0, cos(pitch), -sin(pitch),
				   0, sin(pitch), cos(pitch));

	mat3 Ry = mat3(cos(yaw), 0, sin(yaw),
				   0, 1, 0,
				   -sin(yaw), 0, cos(yaw));

	mat3 Rz = mat3(cos(roll), -sin(roll), 0,
				   sin(roll), cos(roll), 0,
				   0, 0, 1);

	mat3 R = Rx * Ry * Rz;
	vec3 rotated = R * relative;

	// Perspective projection
	p.x = int(focalLength * (rotated.x / rotated.z) + SCREEN_WIDTH / 2);
	p.y = int(focalLength * (rotated.y / rotated.z) + SCREEN_HEIGHT / 2);
}