//DH2323 skeleton code, Lab3 (SDL2 version)
#include <iostream>
#include <glm/glm.hpp>
#include "SDL2auxiliary.h"
#include "TestModel.h"
#include <algorithm>
#include <limits>

using namespace std;
using glm::vec3;
using glm::vec2;
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
vec3 currentColor;

float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

struct Pixel
{
    int x;
    int y;
    float zinv;
};

// ----------------------------------------------------------------------------
// FUNCTIONS
void Update();
void Draw();
void VertexShader(const vec3& v, Pixel& p);
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result);
void DrawLineSDL(ivec2 a, ivec2 b, vec3 color);
void DrawPolygonEdges(const vector<vec3>& vertices);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void DrawPolygon(const vector<vec3>& vertices);

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

void Update()
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

    // Clear the depth buffer at the start of each frame
    for(int y = 0; y < SCREEN_HEIGHT; ++y)
        for(int x = 0; x < SCREEN_WIDTH; ++x)
            depthBuffer[y][x] = 0;

	for (int i = 0; i < triangles.size(); ++i)
	{
		currentColor = triangles[i].color;
		vector<vec3> vertices(3);
		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;
		DrawPolygon(vertices);
	}
	sdlAux->render();
}

// Task 6.4: Modified VertexShader to handle Pixel instead of ivec2
void VertexShader(const vec3& v, Pixel& p)
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

    // Calculate inverse depth
    p.zinv = 1.0f / rotated.z;

	// Perspective projection
	p.x = int(focalLength * rotated.x * p.zinv + SCREEN_WIDTH / 2);
	p.y = int(focalLength * rotated.y * p.zinv + SCREEN_HEIGHT / 2);
}

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result)
{
	int N = result.size();
    if (N <= 0) return;

    float stepX = float(b.x - a.x) / float(max(N - 1, 1));
    float stepY = float(b.y - a.y) / float(max(N - 1, 1));
    float stepZinv = float(b.zinv - a.zinv) / float(max(N - 1, 1));

    float currentX = a.x;
    float currentY = a.y;
    float currentZinv = a.zinv;

	for (int i = 0; i < N; ++i)
	{
		result[i].x = int(currentX);
		result[i].y = int(currentY);
        result[i].zinv = currentZinv;

		currentX += stepX;
		currentY += stepY;
        currentZinv += stepZinv;
	}
}

void DrawLineSDL(ivec2 a, ivec2 b, vec3 color)
{
	ivec2 delta = glm::abs(a - b);
	int pixels = glm::max(delta.x, delta.y) + 1;
	vector<ivec2> line(pixels);

	// Use the vector<ivec2> version of Interpolate
	int N = line.size();
	vec2 step = vec2(b - a) / float(max(N - 1, 1));
	vec2 current(a);
	for (int i = 0; i < N; ++i)
	{
		line[i] = current;
		current += step;
	}

	for (const auto& point : line)
	{
		if (point.x >= 0 && point.x < SCREEN_WIDTH && point.y >= 0 && point.y < SCREEN_HEIGHT)
		{
			sdlAux->putPixel(point.x, point.y, color);
		}
	}
}

void DrawPolygonEdges(const vector<vec3>& vertices)
{
	int V = vertices.size();
	vector<Pixel> projectedVertices(V);
	for (int i = 0; i < V; ++i)
	{
		VertexShader(vertices[i], projectedVertices[i]);
	}
	for (int i = 0; i < V; ++i)
	{
		int j = (i + 1) % V; // The next vertex
		vec3 color(1, 1, 1);
		DrawLineSDL(ivec2(projectedVertices[i].x, projectedVertices[i].y),
                   ivec2(projectedVertices[j].x, projectedVertices[j].y), color);
	}
}

// Task 6.4: Modified ComputePolygonRows to handle Pixel instead of ivec2
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels)
{
	// 1. Find min and max y-value of the polygon
	int minY = vertexPixels[0].y;
	int maxY = vertexPixels[0].y;
	for (const auto& vp : vertexPixels)
	{
		minY = min(minY, vp.y);
		maxY = max(maxY, vp.y);
	}

	// 2. Resize leftPixels and rightPixels
	int numRows = maxY - minY + 1;
	if (numRows <= 0) return; // Handle degenerate cases
	leftPixels.resize(numRows);
	rightPixels.resize(numRows);

	// 3. Initialize the x-coordinates
	for (int i = 0; i < numRows; ++i)
	{
		leftPixels[i].x = numeric_limits<int>::max();
		leftPixels[i].y = minY + i;
        leftPixels[i].zinv = 0;
		rightPixels[i].x = numeric_limits<int>::min();
		rightPixels[i].y = minY + i;
        rightPixels[i].zinv = 0;
	}

	// 4. Loop through all edges and interpolate x-coordinates and zinv
	int V = vertexPixels.size();
	for (int i = 0; i < V; ++i)
	{
		Pixel a = vertexPixels[i];
		Pixel b = vertexPixels[(i + 1) % V];

		int deltaX = abs(a.x - b.x);
		int deltaY = abs(a.y - b.y);
		int steps = glm::max(deltaX, deltaY) + 1;
		vector<Pixel> line(steps);
		Interpolate(a, b, line);

		for (const auto& point : line)
		{
			if (point.y >= minY && point.y <= maxY)
			{
				int row = point.y - minY;
				if (point.x < leftPixels[row].x)
                {
                    leftPixels[row].x = point.x;
                    leftPixels[row].zinv = point.zinv;
                }
                if (point.x > rightPixels[row].x)
                {
                    rightPixels[row].x = point.x;
                    rightPixels[row].zinv = point.zinv;
                }
			}
		}
	}
}

// Task 6.4: Modified DrawRows to handle Pixel instead of ivec2 and use depth buffer
void DrawRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels)
{
	for (size_t i = 0; i < leftPixels.size(); ++i)
	{
		int y = leftPixels[i].y;
        Pixel left = leftPixels[i];
        Pixel right = rightPixels[i];

        // Calculate number of pixels in this row
        int pixels = right.x - left.x + 1;
        if (pixels <= 0) continue;

        // Create vector of pixels for the row
        vector<Pixel> rowPixels(pixels);
        Interpolate(left, right, rowPixels);

		for (const auto& p : rowPixels)
		{
            int x = p.x;
			if (x >= 0 && x < SCREEN_WIDTH && y >= 0 && y < SCREEN_HEIGHT)
			{
                // Check if this pixel is in front of what's already drawn
                if (p.zinv > depthBuffer[y][x])
                {
                    depthBuffer[y][x] = p.zinv; // Update depth buffer
				    sdlAux->putPixel(x, y, currentColor);
                }
			}
		}
	}
}

void DrawPolygon(const vector<vec3>& vertices)
{
	int V = vertices.size();
	vector<Pixel> vertexPixels(V);
	for (int i = 0; i < V; ++i)
	{
		VertexShader(vertices[i], vertexPixels[i]);
	}
	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	DrawRows(leftPixels, rightPixels);
}

// Original Interpolate function kept for backward compatibility with DrawLineSDL
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result)
{
	int N = result.size();
	vec2 step = vec2(b - a) / float(max(N - 1, 1));
	vec2 current(a);
	for (int i = 0; i < N; ++i)
	{
		result[i] = current;
		current += step;
	}
}