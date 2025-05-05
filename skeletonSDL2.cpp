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
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

	const Uint8* keystate = SDL_GetKeyboardState(NULL);
	if (keystate[SDL_SCANCODE_UP]) {
		cameraPos.z += 0.05f;
	}
	if (keystate[SDL_SCANCODE_DOWN]) {
		cameraPos.z -= 0.05f;
	}
	if (keystate[SDL_SCANCODE_LEFT]) {
		cameraPos.x -= 0.05f;
	}
	if (keystate[SDL_SCANCODE_RIGHT]) {
		cameraPos.x += 0.05f;
	}
	if (keystate[SDL_SCANCODE_W]) {

	}
	if (keystate[SDL_SCANCODE_S]) {

	}
	if (keystate[SDL_SCANCODE_A]) {

	}
	if (keystate[SDL_SCANCODE_D]) {

	}
	if (keystate[SDL_SCANCODE_Q]) {

	}
	if (keystate[SDL_SCANCODE_E]) {

	}
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
	vec3 cameraToVertex = v - cameraPos;
	p.x = static_cast<int>(-focalLength * cameraToVertex.x / cameraToVertex.z + SCREEN_WIDTH / 2.0f);
	p.y = static_cast<int>(-focalLength * cameraToVertex.y / cameraToVertex.z + SCREEN_HEIGHT / 2.0f);
}