// DH2323 skeleton code, Lab3 (SDL2 version)
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <glm/glm.hpp>
#include <glm/gtx/constants.hpp>
#include "SDL2auxiliary.h"
#include "TestModel.h"

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

// Camera and viewing parameters
vec3 cameraPos(0, 0, -3.001);
float focalLength = SCREEN_HEIGHT;
float yaw = 0.0f;
float pitch = 0.0f;
float roll = 0.0f;
float rotationSpeed = 0.002f;
float moveSpeed = 0.005f;
float lightMoveSpeed = 0.005f;

float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

vec3 lightPos(0.0f, -0.5f, -0.7f);
vec3 lightPower = 1.0f * vec3(1.0f, 1.0f, 1.0f);
vec3 indirectLightPowerPerArea = 0.3f * vec3(1.0f, 1.0f, 1.0f);

vec3 currentNormal;
vec3 currentReflectance;

struct Vertex
{
    vec3 position;
};

struct Pixel
{
    int x;
    int y;
    float zinv;
    vec3 pos3d;
    vec3 pos3d_over_z;
};


// ----------------------------------------------------------------------------
// FUNCTIONS Prototypes

void Update();
void Draw();
void PixelShader(const Pixel& p);
void VertexShader(const Vertex& v, Pixel& p);
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result);
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void DrawPolygon(const vector<Vertex>& vertices);
void DrawLineSDL(ivec2 a, ivec2 b, vec3 color);
void DrawPolygonEdges(const vector<Vertex>& vertices);

// ----------------------------------------------------------------------------
// MAIN

int SDL_main(int argc, char* argv[])
{
    LoadTestModel(triangles);
    sdlAux = new SDL2Aux(SCREEN_WIDTH, SCREEN_HEIGHT);
    t = SDL_GetTicks(); // Set start value for timer.

    while (!sdlAux->quitEvent())
    {
        Update();
        Draw();
    }

    sdlAux->saveBMP("screenshot_illumination_dim.bmp"); // Changed filename
    delete sdlAux;
    return 0;
}

// ----------------------------------------------------------------------------
// UPDATE

void Update()
{
    int t2 = SDL_GetTicks();
    float dt = float(t2 - t);
    t = t2;
    // cout << "Render time: " << dt << " ms." << endl;

    const Uint8* keystate = SDL_GetKeyboardState(NULL);
    float camMove = moveSpeed * dt;
    float lightMove = lightMoveSpeed * dt;
    float rot = rotationSpeed * dt;

    // --- Camera Movement ---
    vec3 forward(sin(yaw), 0, cos(yaw));
    vec3 right(cos(yaw), 0, -sin(yaw));
    if (keystate[SDL_SCANCODE_UP])    cameraPos += camMove * forward;
    if (keystate[SDL_SCANCODE_DOWN])  cameraPos -= camMove * forward;
    if (keystate[SDL_SCANCODE_LEFT])  cameraPos -= camMove * right;
    if (keystate[SDL_SCANCODE_RIGHT]) cameraPos += camMove * right;

    if (keystate[SDL_SCANCODE_W]) pitch -= rot;
    if (keystate[SDL_SCANCODE_S]) pitch += rot;
    if (keystate[SDL_SCANCODE_Q]) yaw   -= rot;
    if (keystate[SDL_SCANCODE_E]) yaw   += rot;

    // --- Light Movement ---
    if (keystate[SDL_SCANCODE_A]) lightPos.x -= lightMove;
    if (keystate[SDL_SCANCODE_D]) lightPos.x += lightMove;
    if (keystate[SDL_SCANCODE_HOME])     lightPos.y += lightMove;
    if (keystate[SDL_SCANCODE_END])      lightPos.y -= lightMove;
    if (keystate[SDL_SCANCODE_PAGEUP])   lightPos.z += lightMove;
    if (keystate[SDL_SCANCODE_PAGEDOWN]) lightPos.z -= lightMove;

    // Clamp pitch
    float limit = glm::pi<float>() / 2.0f - 0.01f;
    pitch = std::max(-limit, std::min(limit, pitch));
    // roll = 0;
}

// ----------------------------------------------------------------------------
// DRAW

void Draw()
{
    sdlAux->clearPixels();

    // Clear the depth buffer
    for(int y = 0; y < SCREEN_HEIGHT; ++y) {
        for(int x = 0; x < SCREEN_WIDTH; ++x) {
            depthBuffer[y][x] = 0;
        }
    }

    for (size_t i = 0; i < triangles.size(); ++i)
    {
        // Task 7.7: Set global normal and reflectance
        currentNormal = triangles[i].normal;
        currentReflectance = triangles[i].color;

        // Task 7.9: Use simplified Vertex struct
        vector<Vertex> vertices(3);
        vertices[0].position = triangles[i].v0;
        vertices[1].position = triangles[i].v1;
        vertices[2].position = triangles[i].v2;

        DrawPolygon(vertices);
        // DrawPolygonEdges(vertices);
    }

    sdlAux->render();
}


// ----------------------------------------------------------------------------
// SHADERS

// Task 7.9: Modified VertexShader
void VertexShader(const Vertex& v, Pixel& p)
{
    vec3 worldPos = v.position;
    vec3 relative = worldPos - cameraPos;

    mat3 Rx = mat3(1, 0, 0,
                   0, cos(pitch), -sin(pitch),
                   0, sin(pitch), cos(pitch));
    mat3 Ry = mat3(cos(yaw), 0, sin(yaw),
                   0, 1, 0,
                   -sin(yaw), 0, cos(yaw));
    mat3 Rz = mat3(cos(roll), -sin(roll), 0,
                   sin(roll), cos(roll), 0,
                   0, 0, 1);
    mat3 R = Ry * Rx * Rz;

    vec3 rotated = R * relative;

    if (rotated.z <= 0.001f) {
        p.x = -1; p.y = -1; p.zinv = 0; p.pos3d = vec3(0); p.pos3d_over_z = vec3(0); return;
    }

    p.zinv = 1.0f / rotated.z;
    p.x = static_cast<int>(focalLength * rotated.x * p.zinv + SCREEN_WIDTH / 2.0f);
    p.y = static_cast<int>(focalLength * rotated.y * p.zinv + SCREEN_HEIGHT / 2.0f);
    // Y-flip might be needed:
    // p.y = static_cast<int>(-focalLength * rotated.y * p.zinv + SCREEN_HEIGHT / 2.0f);

    p.pos3d = worldPos;
    p.pos3d_over_z = worldPos * p.zinv;
}

// Task 7.9: Modified PixelShader
void PixelShader( const Pixel& p )
{
    int x = p.x;
    int y = p.y;

    if (x < 0 || x >= SCREEN_WIDTH || y < 0 || y >= SCREEN_HEIGHT) {
        return;
    }

    if( p.zinv > depthBuffer[y][x] )
    {
        depthBuffer[y][x] = p.zinv;

        // --- Per-Pixel Illumination Calculation ---
        vec3 worldPos = p.pos3d;
        vec3 normal = glm::normalize(currentNormal);
        vec3 reflectance = currentReflectance;

        vec3 vector_to_light = lightPos - worldPos;
        float distance_to_light_sq = glm::dot(vector_to_light, vector_to_light);
        vec3 dir_to_light_normalized = glm::normalize(vector_to_light);

        float diffuse_intensity = std::max(0.0f, glm::dot(normal, dir_to_light_normalized));

        float attenuation_denominator = distance_to_light_sq;
        if (attenuation_denominator < 1e-5f) attenuation_denominator = 1e-5f;

        vec3 direct_illumination_component = (lightPower * diffuse_intensity) / attenuation_denominator;

        vec3 final_illumination = reflectance * (direct_illumination_component + indirectLightPowerPerArea);
        final_illumination = glm::clamp(final_illumination, 0.0f, 1.0f);

        sdlAux->putPixel( x, y, final_illumination );
    }
}


// ----------------------------------------------------------------------------
// INTERPOLATION

// Task 7.9: Modified Interpolate for perspective-correct pos3d
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result)
{
    size_t N = result.size();
    if (N == 0) return;
    if (N == 1) { result[0] = a; return; }

    float steps = static_cast<float>(N - 1);
    float stepX = static_cast<float>(b.x - a.x) / steps;
    float stepY = static_cast<float>(b.y - a.y) / steps;
    float stepZinv = (b.zinv - a.zinv) / steps;
    vec3 stepPos3dOverZ = (b.pos3d_over_z - a.pos3d_over_z) / steps;

    float currentX = static_cast<float>(a.x);
    float currentY = static_cast<float>(a.y);
    float currentZinv = a.zinv;
    vec3 currentPos3dOverZ = a.pos3d_over_z;

    for (size_t i = 0; i < N; ++i)
    {
        result[i].x = static_cast<int>(round(currentX));
        result[i].y = static_cast<int>(round(currentY));
        result[i].zinv = currentZinv;
        result[i].pos3d_over_z = currentPos3dOverZ;

        if (abs(currentZinv) > 1e-6f) {
            result[i].pos3d = currentPos3dOverZ / currentZinv;
        } else {
             result[i].pos3d = a.pos3d; // Fallback
        }

        currentX += stepX;
        currentY += stepY;
        currentZinv += stepZinv;
        currentPos3dOverZ += stepPos3dOverZ;
    }
}

// Original Interpolate function kept for backward compatibility with DrawLineSDL
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result)
{
    size_t N = result.size();
    if (N == 0) return;
    if (N == 1) { result[0] = a; return; }

    float steps = static_cast<float>(N - 1);
    vec2 step = vec2(b - a) / steps;
    vec2 current(a);

    for (size_t i = 0; i < N; ++i)
    {
        result[i] = ivec2(round(current.x), round(current.y));
        current += step;
    }
}


// ----------------------------------------------------------------------------
// RASTERIZATION

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels)
{
    // 1. Find min/max y
    int minY = SCREEN_HEIGHT, maxY = -1;
    bool validVertexFound = false;
    for (const auto& vp : vertexPixels) {
        if (vp.x != -1) {
             minY = min(minY, vp.y);
             maxY = max(maxY, vp.y);
             validVertexFound = true;
        }
    }
    if (!validVertexFound || maxY < minY) { leftPixels.clear(); rightPixels.clear(); return; }
    minY = max(0, minY);
    maxY = min(SCREEN_HEIGHT - 1, maxY);

    // 2. Resize vectors
    int numRows = maxY - minY + 1;
    if (numRows <= 0) { leftPixels.clear(); rightPixels.clear(); return; }
    leftPixels.resize(numRows);
    rightPixels.resize(numRows);

    // 3. Initialize row endpoints
    for (int i = 0; i < numRows; ++i)
    {
        int currentY = minY + i;
        leftPixels[i].x = numeric_limits<int>::max();
        leftPixels[i].y = currentY;
        leftPixels[i].zinv = 0.0f;
        leftPixels[i].pos3d = vec3(0.0f);
        leftPixels[i].pos3d_over_z = vec3(0.0f);

        rightPixels[i].x = numeric_limits<int>::min();
        rightPixels[i].y = currentY;
        rightPixels[i].zinv = 0.0f;
        rightPixels[i].pos3d = vec3(0.0f);
        rightPixels[i].pos3d_over_z = vec3(0.0f);
    }

    // 4. Interpolate along edges
    int V = vertexPixels.size();
    for (int i = 0; i < V; ++i)
    {
        Pixel a = vertexPixels[i];
        Pixel b = vertexPixels[(i + 1) % V];

        if(a.x == -1 || b.x == -1) continue;

        int dy = abs(a.y - b.y);
        int dx = abs(a.x - b.x);
        int steps = max(dx, dy) + 1;
        steps = max(steps, 1);

        vector<Pixel> edgePixels(steps);
        Interpolate(a, b, edgePixels);

        for (const auto& p_edge : edgePixels)
        {
            if (p_edge.y >= minY && p_edge.y <= maxY)
            {
                int row_index = p_edge.y - minY;
                if (p_edge.x < leftPixels[row_index].x)
                {
                    leftPixels[row_index] = p_edge;
                }
                if (p_edge.x > rightPixels[row_index].x)
                {
                    rightPixels[row_index] = p_edge;
                }
            }
        }
    }
}


void DrawRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels)
{
    if (leftPixels.size() != rightPixels.size() || leftPixels.empty()) return;

    for (size_t i = 0; i < leftPixels.size(); ++i)
    {
        Pixel left = leftPixels[i];
        Pixel right = rightPixels[i];

        if (left.x > right.x) continue;

        int rowWidth = right.x - left.x + 1;
        if (rowWidth <= 0) continue;

        vector<Pixel> rowPixels(rowWidth);
        Interpolate(left, right, rowPixels);

        for (const auto& p_row : rowPixels)
        {
             PixelShader(p_row);
        }
    }
}


void DrawPolygon(const vector<Vertex>& vertices)
{
    int V = vertices.size();
    if (V < 3) return;

    vector<Pixel> vertexPixels(V);
    bool anyVertexValid = false;
    for (int i = 0; i < V; ++i)
    {
        VertexShader(vertices[i], vertexPixels[i]);
        if(vertexPixels[i].x != -1) anyVertexValid = true;
    }

    if (!anyVertexValid) return;

    vector<Pixel> leftPixels;
    vector<Pixel> rightPixels;
    ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
    DrawRows(leftPixels, rightPixels);
}


// ----------------------------------------------------------------------------
// UTILITY FUNCTIONS

void DrawLineSDL(ivec2 a, ivec2 b, vec3 color)
{
    ivec2 delta = glm::abs(a - b);
    int pixels = glm::max(delta.x, delta.y) + 1;
    vector<ivec2> line(pixels);
    Interpolate(a, b, line);
    for (const auto& point : line) {
        if (point.x >= 0 && point.x < SCREEN_WIDTH && point.y >= 0 && point.y < SCREEN_HEIGHT) {
            sdlAux->putPixel(point.x, point.y, color);
        }
    }
}

void DrawPolygonEdges(const vector<Vertex>& vertices)
{
    int V = vertices.size();
    if (V < 2) return;
    vector<Pixel> projectedVertices(V);
    for (int i = 0; i < V; ++i) {
        VertexShader(vertices[i], projectedVertices[i]);
    }
    vec3 white(1.0f, 1.0f, 1.0f);
    for (int i = 0; i < V; ++i) {
        Pixel p1 = projectedVertices[i];
        Pixel p2 = projectedVertices[(i + 1) % V];
        if (p1.x != -1 && p2.x != -1) {
             DrawLineSDL(ivec2(p1.x, p1.y), ivec2(p2.x, p2.y), white);
        }
    }
}
